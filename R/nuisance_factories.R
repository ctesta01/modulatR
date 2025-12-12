# Nuisance Factories for LMTP ----------------------------------------------------

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
make_predictor <- function(
    data,
    formula,
    spec,
    outcome_type = c("continuous", "binary", "density"),
    extra_learner_args = NULL) {
  outcome_type <- match.arg(outcome_type)

  # Case A: Super Learner spec (a list of learners)
  if (is.list(spec) && !is.function(spec)) {
    fit <- nadir::super_learner(
      data = data, formula = formula, learners = spec,
      outcome_type = outcome_type,
      extra_learner_args = extra_learner_args
    )
    return(list(
      kind = "sl",
      fit = fit,
      predict = function(newdata) as.numeric(fit$predict(newdata))
    ))
  }

  # Case B: a specific {nadir} learner constructor is passed
  if ("sl_lnr_name" %in% names(attributes(spec))) {
    predictfun <- do.call(
      what = spec,
      args = c(list(formula = formula, data = data), extra_learner_args))
    return(list(predict = predictfun))
  }

  # Case C: fixed model function, e.g., stats::glm or mgcv::gam
  if (is.function(spec)) {
    fit <- spec(formula = formula, data = data)
    pfun <- function(newdata) {
      as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
    }
    return(list(kind = "fixed", fit = fit, predict = pfun))
  }
  stop("Unsupported model spec: pass a SL library (list) or a fitting function (e.g., glm).")
}

#' Nuisance estimation factory for LMTP
#'
#' @description
#' `LMTPNuisanceFactory` encapsulates the estimation of treatment
#' mechanisms `g_t(A_t | H_t)` and density ratios
#' `r_t(A_t, H_t) = g_t^d(A_t | H_t) / g_t(A_t | H_t)` under an LMTP.
#' It also provides a "K-provider" for computing the cumulative weights
#' `omega_t = prod_{s=1}^t r_s(A_s, H_s)` and a `q_trainer()` used for
#' sequential outcome regression.
#'
#' Two strategies for estimating `r_t` are supported:
#'
#' * `g_mode = "density"` (default): estimate `g_t` via regression and
#'   construct `r_t` via the density ratio, as in Díaz et al. (2023).\
#' * `g_mode = "ratio_classification"`: use the "classification trick"
#'   of Díaz et al. (2023), fitting a classifier to distinguish observed
#'   vs. shifted treatment-history pairs and recovering `r_t` from the
#'   estimated odds.
#'
#' @section References:
#'   Díaz et al. (2023), *Nonparametric Causal Effects Based on
#'   Longitudinal Modified Treatment Policies*.\
#'   Wei et al. (2023), *Efficient Targeted Learning of Heterogeneous
#'   Treatment Effects for Multiple Subgroups*.
#'
#' @export
LMTPNuisanceFactory <- R6::R6Class(
  "LMTPNuisanceFactory",
  public = list(
    # ---- configuration for g ----
    #' @field learners_g List of learner specifications for `g_t` or the
    #'   classifier in the ratio-classification mode.
    learners_g = NULL,
    #' @field learners_r Optional list of learner specifications dedicated
    #'   to the ratio-classification estimator. If `NULL`, `learners_g`
    #'   will be reused.
    learners_r = NULL,
    #' @field fml_g Optional list or single formula for `g_t` models.
    fml_g = NULL,
    #' @field learners_g_extra_args Optional list of extra arguments per
    #'   time point for `learners_g`.
    learners_g_extra_args = NULL,
    #' @field A_type Either "continuous" or "discrete".
    A_type = NULL,
    #' @field policy_seq `LMTPPolicySequence` describing the LMTP.
    policy_seq = NULL,
    #' @field repeat_fmls_lnrs_args Logical: recycle formulas and learner
    #'   specs across time points.
    repeat_fmls_lnrs_args = NULL,
    #' @field g_mode Character: either "density" or "ratio_classification".
    g_mode = NULL,

    # ---- trained artifacts ----
    #' @field g_list List of functions `(A,H) -> g_t(A|H)`; may be `NULL`
    #'   in ratio-classification mode.
    g_list = NULL,
    #' @field r_list List of functions `(A,H) -> r_t(A,H)`.
    r_list = NULL,

    #' @description Create a new LMTP nuisance factory.
    #' @param learners_g Learner spec(s) for the density or classifier.
    #' @param policy_seq `LMTPPolicySequence` describing the LMTP.
    #' @param A_type Either "continuous" or "discrete".
    #' @param fml_g Optional list/single formula for `g_t`.
    #' @param learners_g_extra_args Optional list/single list of extra
    #'   learner args.
    #' @param repeat_fmls_lnrs_args Logical; if `TRUE`, replicate learners
    #'   and formulas over time.
    #' @param g_mode One of "density" (default) or "ratio_classification".
    #' @param learners_r Optional learner spec(s) for the
    #'   ratio-classification estimator.
    initialize = function(learners_g,
                          policy_seq,
                          A_type = c("continuous", "discrete"),
                          fml_g = NULL,
                          learners_g_extra_args = NULL,
                          repeat_fmls_lnrs_args = TRUE,
                          g_mode = c("density", "ratio_classification"),
                          learners_r = NULL) {
      self$learners_g <- learners_g
      self$policy_seq <- policy_seq
      self$A_type <- match.arg(A_type)
      self$fml_g <- fml_g
      self$learners_g_extra_args <- learners_g_extra_args
      self$repeat_fmls_lnrs_args <- repeat_fmls_lnrs_args
      self$g_mode <- match.arg(g_mode)
      self$learners_r <- learners_r
    },

    #' @description Train `g_t` / `r_t` across all time points.
    #' @param ds An `LMTP_Data_Struct` instance.
    #' @return Invisibly returns `self` with `g_list` and `r_list` populated.
    train = function(ds) {
      tau <- ds$tau()
      learners_g <- .rep_if_needed(self$learners_g, tau, self$repeat_fmls_lnrs_args)
      learners_g_extra_args <- .rep_if_needed(self$learners_g_extra_args, tau, self$repeat_fmls_lnrs_args)
      fml_g <- if (is.null(self$fml_g)) vector("list", tau) else
        .rep_if_needed(self$fml_g, tau, self$repeat_fmls_lnrs_args)

      g_list <- vector("list", tau)
      r_list <- vector("list", tau)

      for (t in seq_len(tau)) {
        res <- local({
          t_i <- t
          AH_i <- ds$AH(t_i)
          H_i  <- ds$H(t_i)
          Aname_i <- ds$A_cols[[t_i]]

          if (is.null(fml_g[[t_i]])) {
            fml_g[[t_i]] <- stats::as.formula(
              paste0(Aname_i, " ~ ", paste0(colnames(H_i), collapse = " + "))
            )
          }

          if (identical(self$g_mode, "density")) {
            # Density-based g_t and ratio r_t = g_t^d / g_t
            sl_g_i <- make_predictor(
              data = AH_i, formula = fml_g[[t_i]],
              spec = learners_g[[t_i]],
              outcome_type = "density",
              extra_learner_args = learners_g_extra_args[[t_i]]
            )

            g_pred_fun_i <- function(A_vec, H_df) {
              nd <- dplyr::mutate(H_df, "{Aname_i}" := A_vec)
              as.numeric(sl_g_i$predict(nd))
            }

            if (self$A_type == "continuous") {
              r_fun_i <- function(new_A, new_H) {
                g_obs <- g_pred_fun_i(new_A, new_H)
                gd_obs <- self$policy_seq$gd_from_density_t(
                  t_i, new_A, new_H, density_fun = g_pred_fun_i
                )
                pmax(gd_obs, 1e-12) / pmax(g_obs, 1e-12)
              }
            } else {
              # binary / discrete A
              r_fun_i <- function(new_A, new_H) {
                n <- nrow(new_H)
                A0 <- rep(0, n); A1 <- rep(1, n)
                dA0 <- self$policy_seq$apply_policy_t(t_i, A0, new_H)
                dA1 <- self$policy_seq$apply_policy_t(t_i, A1, new_H)
                g1 <- g_pred_fun_i(A1, new_H); g0 <- g_pred_fun_i(A0, new_H)
                gd1 <- (dA0 == 1) * g0 + (dA1 == 1) * g1
                gd0 <- (dA0 == 0) * g0 + (dA1 == 0) * g1
                gobs  <- ifelse(new_A == 1, g1, g0)
                gdobs <- ifelse(new_A == 1, gd1, gd0)
                pmax(gdobs, 1e-12) / pmax(gobs, 1e-12)
              }
            }

            list(g = g_pred_fun_i, r = r_fun_i)
          } else {
            # ratio_classification: classification trick for r_t
            A_obs <- ds$A(t_i)
            A_shift <- self$policy_seq$apply_policy_t(t_i, A_obs, H_i)

            df0 <- dplyr::mutate(AH_i, .lmpt_shift = 0L)
            df1 <- dplyr::mutate(dplyr::mutate(H_i, "{Aname_i}" := A_shift), .lmpt_shift = 1L)
            data_cls <- dplyr::bind_rows(df0, df1)

            fml_cls <- stats::as.formula(
              paste0(".lmpt_shift ~ ", paste0(c(colnames(H_i), Aname_i), collapse = " + "))
            )

            learners_r <- if (!is.null(self$learners_r)) .rep_if_needed(self$learners_r, tau, self$repeat_fmls_lnrs_args) else learners_g

            sl_r_i <- make_predictor(
              data = data_cls,
              formula = fml_cls,
              spec = learners_r[[t_i]],
              outcome_type = "binary",
              extra_learner_args = learners_g_extra_args[[t_i]]
            )

            r_fun_i <- function(new_A, new_H) {
              nd <- dplyr::mutate(new_H, "{Aname_i}" := new_A)
              p1 <- as.numeric(sl_r_i$predict(nd))
              p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)
              p1 / (1 - p1)  # odds = g^d / g
            }

            g_dummy <- function(...) {
              stop("g_t is not available in g_mode = 'ratio_classification'; only r_t is estimated.")
            }

            list(g = g_dummy, r = r_fun_i)
          }
        })

        g_list[[t]] <- res$g
        r_list[[t]] <- res$r
      }

      self$g_list <- g_list
      self$r_list <- r_list
      invisible(self)
    },

    #' @description Build a provider for cumulative LMTP weights.
    #' @param ds `LMTP_Data_Struct` instance.
    #' @return A list with two functions:
    #'   * `K_obs(t)` returns the length-n vector of `omega_t`.
    #'   * `K_obs_all` is the list of all `omega_t` vectors.
    k_provider = function(ds) {
      if (is.null(self$r_list)) stop("Call $train(ds) before requesting K-provider.")
      tau <- ds$tau()
      n <- nrow(ds$data)
      K_obs_t <- vector("list", length = tau)

      for (t in seq_len(tau)) {
        K <- rep(1, n)
        for (k in seq_len(t)) {
          K <- K * self$r_list[[k]](ds$A(k), ds$H(k))
        }
        K_obs_t[[t]] <- K
      }

      list(
        K_obs = function(t) K_obs_t[[t]],
        K_obs_all = K_obs_t
      )
    },

    #' @description Build a `Q`-trainer for sequential outcome regression.
    #' @param learners_Q List of learner specifications for `m_t`.
    #' @param fml_Q Optional formula or list of formulas for the outcome
    #'   regressions.
    #' @param repeat_fmls Logical; if `TRUE`, reuse the same formula across
    #'   time points when only one is supplied.
    #' @return A function `train_Q_t(ds, t, pseudo_outcome_vec)` that fits
    #'   `m_t` and returns a prediction function `(A, H) -> m_t(A, H)`.
    q_trainer = function(learners_Q, fml_Q = NULL, repeat_fmls = TRUE) {
      function(ds, t, pseudo_outcome_vec) {
        Aname <- ds$A_cols[[t]]; Ht <- ds$H(t); AHt <- ds$AH(t)
        tau <- ds$tau()

        fmls <- if (is.null(fml_Q)) vector("list", tau) else {
          if (repeat_fmls) rep(list(fml_Q), tau) else fml_Q
        }
        if (is.null(fmls[[t]])) {
          rhs <- paste0(c(colnames(Ht), Aname), collapse = " + ")
          fmls[[t]] <- stats::as.formula(paste0("M_next ~ ", rhs))
        }

        data_fit <- dplyr::bind_cols(AHt, M_next = pseudo_outcome_vec)
        sl_Q <- make_predictor(
          data = data_fit, formula = fmls[[t]],
          spec = learners_Q[[t]], outcome_type = "continuous"
        )

        function(A_vec, H_df) {
          nd <- dplyr::bind_cols(H_df, !!rlang::sym(Aname) := A_vec, M_next = NA_real_)
          as.numeric(sl_Q$predict(nd))
        }
      }
    }
  )
)
