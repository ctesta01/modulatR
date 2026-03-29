
# Nuisance Factories for LMTP ------------------------------------------------

#' LMTP nuisance estimation factory
#'
#' @description
#' `LMTPNuisanceFactory` encapsulates estimation of treatment mechanisms
#' `g_t(A_t | H_t)` and density ratios
#' `r_t(A_t, H_t) = g_t^d(A_t | H_t) / g_t(A_t | H_t)` under an LMTP.
#'
#' It also provides:
#' \itemize{
#'   \item a provider for cumulative weights
#'   `omega_t = \prod_{s=1}^t r_s(A_s, H_s)`,
#'   \item a `q_trainer()` helper for sequential outcome regression.
#' }
#'
#' Two strategies for estimating `r_t` are supported:
#' \itemize{
#'   \item `g_mode = "density"`: estimate `g_t` and form
#'   `r_t = g_t^d / g_t`;
#'   \item `g_mode = "ratio_classification"`: estimate `r_t`
#'   directly via the classification trick.
#' }
#'
#' @export
#' @examples
#' #' @examples
#' @examples
#' ## -------------------------------------------------
#' ## Example: LMTPNuisanceFactory in density mode
#' ## -------------------------------------------------
#' set.seed(1)
#' n <- 200
#'
#' df <- data.frame(
#'   W1 = rnorm(n),
#'   L1 = rnorm(n)
#' )
#' df$A1 <- rbinom(n, 1, plogis(0.6 * df$W1 - 0.8 * df$L1))
#' df$L2 <- rnorm(n, mean = 0.3 * df$A1 + 0.2 * df$W1)
#' df$A2 <- rbinom(n, 1, plogis(-0.4 * df$W1 + 0.7 * df$L2))
#' df$Y  <- rnorm(n, mean = 0.5 * df$A1 + 0.8 * df$A2 + 0.3 * df$L2)
#' ds <- LMTPData$new(
#'   data = df,
#'   A_cols = c("A1", "A2"),
#'   L_cols = list("L1", "L2"),
#'   W_cols = "W1",
#'   Y_col = "Y"
#' )
#'
#' pol_flip <- mtp_discrete(
#'   map_fun = function(A, H) 1 - A,
#'   support = c(0, 1)
#' )
#' policy_seq <- repeat_policy_over_time(pol_flip, tau = ds$tau())
#'
#' nf <- LMTPNuisanceFactory$new(
#'   learners_g = stats::glm,
#'   policy_seq = policy_seq,
#'   A_type = "discrete",
#'   g_mode = "density",
#'   learners_g_extra_args = list(family = stats::binomial())
#' )
#'
#' nf$train(ds)
#' nf
#'
#' head(nf$g_list[[1]](ds$A(1), ds$H(1)))
#' head(nf$r_list[[1]](ds$A(1), ds$H(1)))
#'
#' kp <- nf$k_provider(ds)
#' head(kp$K_obs(1))
#' head(kp$K_obs(2))
#'
#' ## -------------------------------------------------
#' ## Example: q_trainer() for sequential regression
#' ## -------------------------------------------------
#' q_fit_factory <- nf$q_trainer(
#'   learners_Q = stats::glm
#' )
#'
#' pseudo_outcome <- ds$Y()
#' Q2 <- q_fit_factory(ds, t = 2, pseudo_outcome_vec = pseudo_outcome)
#'
#' head(Q2(ds$A(2), ds$H(2)))
LMTPNuisanceFactory <- R6::R6Class(
  "LMTPNuisanceFactory",
  public = list(
    # ---- configuration ----
    learners_g = NULL,
    learners_r = NULL,
    fml_g = NULL,
    learners_g_extra_args = NULL,
    A_type = NULL,
    policy_seq = NULL,
    repeat_fmls_lnrs_args = NULL,
    g_mode = NULL,

    cross_fit = FALSE,
    V = 5,
    fold_id = NULL,
    seed = NULL,

    g_fits = NULL,
    r_fits = NULL,
    g_oof = NULL,
    r_oof = NULL,
    K_oof = NULL,
    # ---- trained artifacts ----
    g_list = NULL,
    r_list = NULL,

    #' @description Create a new nuisance factory.
    #' @param learners_g Learner specification(s) for `g_t` or the classifier.
    #' @param policy_seq An `LMTPPolicySequence`.
    #' @param A_type Either `"continuous"` or `"discrete"`.
    #' @param fml_g Optional formula or list of formulas for `g_t`.
    #' @param learners_g_extra_args Optional extra learner args.
    #' @param repeat_fmls_lnrs_args Logical; recycle length-1 inputs over time.
    #' @param g_mode Either `"density"` or `"ratio_classification"`.
    #' @param learners_r Optional learner specification(s) for ratio classification.
    initialize = function(learners_g,
                          policy_seq,
                          A_type = c("continuous", "discrete"),
                          fml_g = NULL,
                          learners_g_extra_args = NULL,
                          repeat_fmls_lnrs_args = TRUE,
                          g_mode = c("density", "ratio_classification"),
                          learners_r = NULL,
                          cross_fit = FALSE,
                          V = 5,
                          fold_id = NULL,
                          seed = NULL) {
      self$learners_g <- learners_g
      self$policy_seq <- policy_seq
      self$A_type <- match.arg(A_type)
      self$fml_g <- fml_g
      self$learners_g_extra_args <- learners_g_extra_args
      self$repeat_fmls_lnrs_args <- repeat_fmls_lnrs_args
      self$g_mode <- match.arg(g_mode)
      self$learners_r <- learners_r

      self$cross_fit <- cross_fit
      self$V <- V
      self$fold_id <- fold_id
      self$seed <- seed

      private$validate()
      invisible(self)
    },

    #' @description Train `g_t` / `r_t` across time points.
    #' @param ds An `LMTPData` instance.
    #' @return Invisibly returns `self`.
    train = function(ds) {
      if (!inherits(ds, "LMTPData")) {
        stop("`ds` must inherit from `LMTPData`.")
      }

      self$policy_seq$validate_against_data(ds)

      tau <- ds$tau()
      n <- ds$n

      learners_g <- .rep_if_needed(
        self$learners_g, tau, self$repeat_fmls_lnrs_args
      )
      learners_g_extra_args <- .rep_if_needed(
        self$learners_g_extra_args, tau, self$repeat_fmls_lnrs_args
      )

      fml_g <- if (is.null(self$fml_g)) {
        vector("list", tau)
      } else {
        .rep_if_needed(self$fml_g, tau, self$repeat_fmls_lnrs_args)
      }

      learners_r <- if (!is.null(self$learners_r)) {
        .rep_if_needed(self$learners_r, tau, self$repeat_fmls_lnrs_args)
      } else {
        learners_g
      }

      g_list <- vector("list", tau)
      r_list <- vector("list", tau)
      g_fits <- vector("list", tau)
      r_fits <- vector("list", tau)

      # ------------------------------------------------------------------
      # No cross-fitting
      # ------------------------------------------------------------------
      if (!isTRUE(self$cross_fit)) {
        for (t in seq_len(tau)) {
          AH_t <- ds$AH(t)
          H_t <- ds$H(t)
          Aname_t <- ds$A_cols[[t]]

          if (is.null(fml_g[[t]])) {
            rhs <- setdiff(colnames(H_t), ds$id_col)
            fml_g[[t]] <- stats::as.formula(
              paste0(Aname_t, " ~ ", paste(rhs, collapse = " + "))
            )
          }

          if (identical(self$g_mode, "density")) {
            pred_g_t <- make_predictor(
              data = AH_t,
              formula = fml_g[[t]],
              spec = learners_g[[t]],
              outcome_type = if (self$A_type == "continuous") "density" else "binary",
              extra_learner_args = learners_g_extra_args[[t]]
            )

            g_fits[[t]] <- pred_g_t

            if (identical(self$A_type, "continuous")) {
              g_pred_fun_t <- local({
                pred_local <- pred_g_t
                Aname_local <- Aname_t

                function(A_vec, H_df) {
                  nd <- H_df
                  nd[[Aname_local]] <- A_vec
                  as.numeric(pred_local$predict(nd))
                }
              })
            } else {
              g_pred_fun_t <- local({
                pred_local <- pred_g_t
                Aname_local <- Aname_t

                function(A_vec, H_df) {
                  if (!all(A_vec %in% c(0, 1))) {
                    stop("Discrete density mode currently assumes binary treatment coded 0/1.")
                  }

                  nd1 <- H_df
                  nd1[[Aname_local]] <- 1
                  p1 <- as.numeric(pred_local$predict(nd1))
                  p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)

                  ifelse(A_vec == 1, p1, 1 - p1)
                }
              })
            }

            if (identical(self$A_type, "continuous")) {
              r_fun_t <- local({
                g_local <- g_pred_fun_t
                t_local <- t

                function(new_A, new_H) {
                  g_obs <- g_local(new_A, new_H)
                  gd_obs <- self$policy_seq$gd_t(
                    t = t_local,
                    A_vec = new_A,
                    H_df = new_H,
                    density_fun = g_local
                  )
                  pmax(gd_obs, 1e-12) / pmax(g_obs, 1e-12)
                }
              })
            } else {
              r_fun_t <- local({
                g_local <- g_pred_fun_t
                t_local <- t

                function(new_A, new_H) {
                  if (!all(new_A %in% c(0, 1))) {
                    stop("Discrete density mode currently assumes binary treatment coded 0/1.")
                  }

                  n_local <- nrow(new_H)
                  A0 <- rep(0, n_local)
                  A1 <- rep(1, n_local)

                  dA0 <- self$policy_seq$apply_t(t_local, A0, new_H)
                  dA1 <- self$policy_seq$apply_t(t_local, A1, new_H)

                  g1 <- g_local(A1, new_H)
                  g0 <- g_local(A0, new_H)

                  gd1 <- (dA0 == 1) * g0 + (dA1 == 1) * g1
                  gd0 <- (dA0 == 0) * g0 + (dA1 == 0) * g1

                  gobs <- ifelse(new_A == 1, g1, g0)
                  gdobs <- ifelse(new_A == 1, gd1, gd0)

                  pmax(gdobs, 1e-12) / pmax(gobs, 1e-12)
                }
              })
            }

            g_list[[t]] <- g_pred_fun_t
            r_list[[t]] <- r_fun_t
          } else {
            A_obs <- ds$A(t)
            A_shift <- self$policy_seq$apply_t(t, A_obs, H_t)

            df0 <- AH_t
            df0$.lmpt_shift <- 0L

            df1 <- H_t
            df1[[Aname_t]] <- A_shift
            df1$.lmpt_shift <- 1L

            data_cls <- rbind(df0, df1)

            rhs_cls <- setdiff(colnames(df1), c(".lmpt_shift", ds$id_col))
            fml_cls <- stats::as.formula(
              paste0(".lmpt_shift ~ ", paste(rhs_cls, collapse = " + "))
            )

            pred_r_t <- make_predictor(
              data = data_cls,
              formula = fml_cls,
              spec = learners_r[[t]],
              outcome_type = "binary",
              extra_learner_args = learners_g_extra_args[[t]]
            )

            r_fits[[t]] <- pred_r_t

            r_fun_t <- local({
              pred_local <- pred_r_t
              Aname_local <- Aname_t

              function(new_A, new_H) {
                nd <- new_H
                nd[[Aname_local]] <- new_A
                p1 <- as.numeric(pred_local$predict(nd))
                p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)
                p1 / (1 - p1)
              }
            })

            g_dummy <- function(...) {
              stop("g_t is not available when `g_mode = 'ratio_classification'`; only r_t is estimated.")
            }

            g_list[[t]] <- g_dummy
            r_list[[t]] <- r_fun_t
          }
        }

        self$g_list <- g_list
        self$r_list <- r_list
        self$g_fits <- g_fits
        self$r_fits <- r_fits
        self$g_oof <- NULL
        self$r_oof <- NULL
        self$K_oof <- NULL

        return(invisible(self))
      }

      # ------------------------------------------------------------------
      # Cross-fitting
      # ------------------------------------------------------------------
      fold_id <- .make_folds(
        n = n,
        V = self$V,
        fold_id = self$fold_id,
        seed = self$seed
      )
      self$fold_id <- fold_id

      g_oof <- vector("list", tau)
      r_oof <- vector("list", tau)

      for (t in seq_len(tau)) {
        g_oof[[t]] <- rep(NA_real_, n)
        r_oof[[t]] <- rep(NA_real_, n)

        g_fits[[t]] <- vector("list", length(unique(fold_id)))
        r_fits[[t]] <- vector("list", length(unique(fold_id)))

        Aname_t <- ds$A_cols[[t]]
        H_t_full <- ds$H(t)

        if (is.null(fml_g[[t]])) {
          rhs <- setdiff(colnames(H_t_full), ds$id_col)
          fml_g[[t]] <- stats::as.formula(
            paste0(Aname_t, " ~ ", paste(rhs, collapse = " + "))
          )
        }

        for (v in sort(unique(fold_id))) {
          train_idx <- which(fold_id != v)
          valid_idx <- which(fold_id == v)

          ds_tr <- .subset_ds(ds, train_idx)

          H_valid <- ds$H(t)[valid_idx, , drop = FALSE]
          A_valid <- ds$A(t)[valid_idx]

          if (identical(self$g_mode, "density")) {
            pred_g_t <- make_predictor(
              data = ds_tr$AH(t),
              formula = fml_g[[t]],
              spec = learners_g[[t]],
              outcome_type = if (self$A_type == "continuous") "density" else "binary",
              extra_learner_args = learners_g_extra_args[[t]]
            )

            g_fits[[t]][[v]] <- pred_g_t

            if (identical(self$A_type, "continuous")) {
              g_valid_fun <- local({
                pred_local <- pred_g_t
                Aname_local <- Aname_t

                function(A_vec, H_df) {
                  nd <- H_df
                  nd[[Aname_local]] <- A_vec
                  as.numeric(pred_local$predict(nd))
                }
              })

              g_oof[[t]][valid_idx] <- g_valid_fun(A_valid, H_valid)

              gd_valid <- self$policy_seq$gd_t(
                t = t,
                A_vec = A_valid,
                H_df = H_valid,
                density_fun = g_valid_fun
              )

              r_oof[[t]][valid_idx] <- pmax(gd_valid, 1e-12) / pmax(g_oof[[t]][valid_idx], 1e-12)
            } else {
              g_valid_fun <- local({
                pred_local <- pred_g_t
                Aname_local <- Aname_t

                function(A_vec, H_df) {
                  if (!all(A_vec %in% c(0, 1))) {
                    stop("Discrete density mode currently assumes binary treatment coded 0/1.")
                  }

                  nd1 <- H_df
                  nd1[[Aname_local]] <- 1
                  p1 <- as.numeric(pred_local$predict(nd1))
                  p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)

                  ifelse(A_vec == 1, p1, 1 - p1)
                }
              })

              g_oof[[t]][valid_idx] <- g_valid_fun(A_valid, H_valid)

              n_valid <- length(valid_idx)
              A0 <- rep(0, n_valid)
              A1 <- rep(1, n_valid)

              dA0 <- self$policy_seq$apply_t(t, A0, H_valid)
              dA1 <- self$policy_seq$apply_t(t, A1, H_valid)

              g1 <- g_valid_fun(A1, H_valid)
              g0 <- g_valid_fun(A0, H_valid)

              gd1 <- (dA0 == 1) * g0 + (dA1 == 1) * g1
              gd0 <- (dA0 == 0) * g0 + (dA1 == 0) * g1

              gobs <- ifelse(A_valid == 1, g1, g0)
              gdobs <- ifelse(A_valid == 1, gd1, gd0)

              r_oof[[t]][valid_idx] <- pmax(gdobs, 1e-12) / pmax(gobs, 1e-12)
            }
          } else {
            A_tr <- ds_tr$A(t)
            H_tr <- ds_tr$H(t)
            A_shift_tr <- self$policy_seq$apply_t(t, A_tr, H_tr)

            df0 <- ds_tr$AH(t)
            df0$.lmpt_shift <- 0L

            df1 <- H_tr
            df1[[Aname_t]] <- A_shift_tr
            df1$.lmpt_shift <- 1L

            data_cls <- rbind(df0, df1)

            rhs_cls <- setdiff(colnames(df1), c(".lmpt_shift", ds$id_col))
            fml_cls <- stats::as.formula(
              paste0(".lmpt_shift ~ ", paste(rhs_cls, collapse = " + "))
            )

            pred_r_t <- make_predictor(
              data = data_cls,
              formula = fml_cls,
              spec = learners_r[[t]],
              outcome_type = "binary",
              extra_learner_args = learners_g_extra_args[[t]]
            )

            r_fits[[t]][[v]] <- pred_r_t

            nd <- H_valid
            nd[[Aname_t]] <- A_valid
            p1 <- as.numeric(pred_r_t$predict(nd))
            p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)

            r_oof[[t]][valid_idx] <- p1 / (1 - p1)
          }
        }
      }

      K_oof <- vector("list", tau)
      for (t in seq_len(tau)) {
        K <- rep(1, n)
        for (k in seq_len(t)) {
          K <- K * r_oof[[k]]
        }
        K_oof[[t]] <- K
      }

      self$g_list <- NULL
      self$r_list <- NULL
      self$g_fits <- g_fits
      self$r_fits <- r_fits
      self$g_oof <- g_oof
      self$r_oof <- r_oof
      self$K_oof <- K_oof

      return(invisible(self))
    },

    #' @description Build a provider for cumulative weights.
    #' @param ds An `LMTPData` instance.
    #' @return A list with `K_obs(t)` and `K_obs_all`.
    k_provider = function(ds) {
      if (!inherits(ds, "LMTPData")) {
        stop("`ds` must inherit from `LMTPData`.")
      }

      tau <- ds$tau()

      if (!is.null(self$K_oof)) {
        K_obs_t <- self$K_oof

        return(list(
          K_obs = function(t) {
            if (!is.numeric(t) || length(t) != 1L || t < 1L || t > tau) {
              stop("`t` must be an integer between 1 and tau.")
            }
            K_obs_t[[t]]
          },
          K_obs_all = K_obs_t
        ))
      }

      if (is.null(self$r_list)) {
        stop("Call `$train(ds)` before requesting the K-provider.")
      }

      n <- ds$n
      K_obs_t <- vector("list", length = tau)

      for (t in seq_len(tau)) {
        K <- rep(1, n)
        for (k in seq_len(t)) {
          K <- K * self$r_list[[k]](ds$A(k), ds$H(k))
        }
        K_obs_t[[t]] <- K
      }

      list(
        K_obs = function(t) {
          if (!is.numeric(t) || length(t) != 1L || t < 1L || t > tau) {
            stop("`t` must be an integer between 1 and tau.")
          }
          K_obs_t[[t]]
        },
        K_obs_all = K_obs_t
      )
    },

    #' @description Build a trainer for sequential outcome regressions.
    #' @param learners_Q Learner specification(s) for `Q_t`.
    #' @param fml_Q Optional formula or list of formulas.
    #' @param repeat_fmls Logical; whether to recycle length-1 formulas.
    #' @param learners_Q_extra_args Optional extra args for `Q` learners.
    #' @return A function `train_Q_t(ds, t, pseudo_outcome_vec)`.
    q_trainer = function(learners_Q,
                         fml_Q = NULL,
                         repeat_fmls = TRUE,
                         learners_Q_extra_args = NULL) {

      force(learners_Q)
      force(fml_Q)
      force(repeat_fmls)
      force(learners_Q_extra_args)

      function(ds, t, pseudo_outcome_vec) {
        if (!inherits(ds, "LMTPData")) {
          stop("`ds` must inherit from `LMTPData`.")
        }

        tau <- ds$tau()

        learners_Q_t <- .rep_if_needed(learners_Q, tau, repeat_fmls)

        fmls_Q_t <- if (is.null(fml_Q)) {
          vector("list", tau)
        } else {
          .rep_if_needed(fml_Q, tau, repeat_fmls)
        }

        learners_Q_extra_args_t <- .rep_if_needed(
          learners_Q_extra_args, tau, repeat_fmls
        )

        Aname <- ds$A_cols[[t]]
        Ht <- ds$H(t)
        AHt <- ds$AH(t)

        if (is.null(fmls_Q_t[[t]])) {
          rhs <- setdiff(colnames(AHt), c("M_next", ds$id_col))
          fmls_Q_t[[t]] <- stats::as.formula(
            paste0("M_next ~ ", paste(rhs, collapse = " + "))
          )
        }

        data_fit <- AHt
        data_fit$M_next <- pseudo_outcome_vec

        pred_Q_t <- make_predictor(
          data = data_fit,
          formula = fmls_Q_t[[t]],
          spec = learners_Q_t[[t]],
          outcome_type = "continuous",
          extra_learner_args = learners_Q_extra_args_t[[t]]
        )

        function(A_vec, H_df) {
          nd <- H_df
          nd[[Aname]] <- A_vec
          as.numeric(pred_Q_t$predict(nd))
        }
      }
    },

    #' @description Print a compact summary.
    #' @return The object invisibly.
    print = function(...) {
      cat("LMTPNuisanceFactory\n")
      cat("  tau: ", self$policy_seq$tau(), "\n", sep = "")
      cat("  A_type: ", self$A_type, "\n", sep = "")
      cat("  g_mode: ", self$g_mode, "\n", sep = "")
      cat("  trained: ", !is.null(self$r_list), "\n", sep = "")
      invisible(self)
    },

    get_fold_id = function() {
      self$fold_id
    }
  ),

  private = list(
    validate = function() {
      if (!inherits(self$policy_seq, "LMTPPolicySequence")) {
        stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
      }

      if (!is.logical(self$cross_fit) || length(self$cross_fit) != 1L) {
        stop("`cross_fit` must be TRUE or FALSE.")
      }

      if (!is.numeric(self$V) || length(self$V) != 1L || self$V < 2) {
        stop("`V` must be a single integer >= 2.")
      }

      if (!is.null(self$fold_id) && (!is.numeric(self$fold_id) || any(self$fold_id < 1))) {
        stop("`fold_id` must be NULL or a positive integer vector.")
      }
    }
  )
)





#' Build a prediction function from a learner specification
#'
#' @description
#' Helper that normalizes a variety of learner specifications into a
#' common interface returning a list with a `$predict(newdata)` method.
#' This is used for both outcome regressions and treatment models.
#'
#' @param data Data frame used for training.
#' @param formula Model formula.
#' @param spec Either a list of learners for `nadir::super_learner`, a
#'   single nadir learner function, or a base R fitting function (e.g.,
#'   [stats::glm]).
#' @param outcome_type One of `"continuous"`, `"binary"`, or `"density"`.
#' @param extra_learner_args Optional list of extra arguments passed to
#'   the learner.
#'
#' @return A list with components `kind`, `fit`, and a function `predict(newdata)`.
#' @export
#' @examples
#' ## -------------------------------------------------
#' ## Example: make_predictor() with a GLM
#' ## -------------------------------------------------
#' set.seed(1)
#' n <- 100
#'
#' df <- data.frame(
#'   W1 = rnorm(n),
#'   W2 = rnorm(n)
#' )
#' df$A <- rbinom(n, 1, plogis(0.5 * df$W1 - 0.25 * df$W2))
#'
#' pred <- make_predictor(
#'   data = df,
#'   formula = A ~ W1 + W2,
#'   spec = stats::glm,
#'   outcome_type = "binary",
#'   extra_learner_args = list(family = stats::binomial())
#' )
#'
#' head(pred$predict(df))
is_nadir_learner <- function(x) {
  is.function(x) && "sl_lnr_name" %in% names(attributes(x))
}


# helpers for LMTPNuisanceFactory -----------------------------------------


make_predictor <- function(
    data,
    formula,
    spec,
    outcome_type = c("continuous", "binary", "density"),
    extra_learner_args = NULL) {
  outcome_type <- match.arg(outcome_type)

  # Case A: Super Learner spec = list of nadir learners
  if (is.list(spec) && !is.function(spec)) {
    if (!all(vapply(spec, is_nadir_learner, logical(1)))) {
      stop(
        "If `spec` is a list, it must be a list of `{nadir}` learners ",
        "for `nadir::super_learner()`. ",
        "For a single model like `stats::glm`, pass the function directly, not inside a list."
      )
    }

    fit <- nadir::super_learner(
      data = data,
      formulas = list(.default = formula),
      learners = spec,
      outcome_type = outcome_type,
      extra_learner_args = extra_learner_args
    )

    return(list(
      kind = "sl",
      fit = fit,
      predict = function(newdata) as.numeric(stats::predict(fit, newdata = newdata))
    ))
  }

  # Case B: a single nadir learner constructor
  if (is_nadir_learner(spec)) {
    predictfun <- do.call(
      what = spec,
      args = c(list(formula = formula, data = data), extra_learner_args)
    )

    return(list(
      kind = "nadir_learner",
      fit = NULL,
      predict = function(newdata) as.numeric(predictfun(newdata))
    ))
  }

  # Case C: fixed model function, e.g. stats::glm
  if (is.function(spec)) {
    fit <- do.call(
      what = spec,
      args = c(list(formula = formula, data = data), extra_learner_args)
    )

    return(list(
      kind = "fixed",
      fit = fit,
      predict = function(newdata) {
        as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
      }
    ))
  }

  stop(
    "Unsupported model spec: pass either a single fitting function, ",
    "a single `{nadir}` learner, or a list of `{nadir}` learners."
  )
}

