#' Internal utility: replicate a list or object over time
#'
#' Helper used throughout the LMTP implementation to expand a single
#' specification (e.g., a learner list or formula) to a list of length
#' `tau`, when `repeat_bool = TRUE`. If `repeat_bool = FALSE`, the input
#' must already be a list of length `tau`.
#'
#' @param x An object or list to be replicated.
#' @param tau Integer number of time points.
#' @param repeat_bool Logical; if `TRUE`, replicate `x` `tau` times; if
#'   `FALSE`, check that `x` is already a list of length `tau`.
#'
#' @return A list of length `tau`.
.rep_if_needed <- function(x, tau, repeat_bool) {

  if (repeat_bool) {
    x <- rep(list(x), tau)
  } else {
    if (length(x) != tau || !is.list(x)) {
      stop("formulas, learner lists, and learner arguments must be length tau if not using repeat_fmls_lnrs_and_args")
    }
  }
  return(x)
}

# Policy Sequence Class -------------------------------------------------------

#' LMTP policy sequence
#'
#' @description
#' `LMTPPolicySequence` stores a list of single-time modified treatment
#' policies (MTPs) and provides convenience methods to apply them over
#' time. Each element of `policies` must implement at least two methods:
#' `apply_policy(A_vec, H_df)` and `gd_from_density(A_vec, H_df, density_fun)`.
#'
#' @section References:
#'   Longitudinal modified treatment policies are described in Diaz et al. (2023).\
#'   See also the `lmtp` package for reference implementations.
#'
#' @export
LMTPPolicySequence <- R6::R6Class(
  "LMTPPolicySequence",
  public = list(
    #' @field policies List of per-time MTP objects (length tau).
    policies = NULL,
    #' @description Create a new LMTP policy sequence.
    #' @param policies List of per-time MTP objects.
    initialize = function(policies) {
      stopifnot(is.list(policies), length(policies) >= 1)
      self$policies <- policies
    },
    #' @description Number of time points (tau).
    #' @return Integer number of time points.
    tau = function() length(self$policies),
    #' @description Apply the t-th policy to a treatment vector.
    #' @param t Time index (1-based).
    #' @param A_vec Numeric or factor vector of observed treatment at time `t`.
    #' @param H_df Data frame with history covariates `H_t`.
    #' @return Vector of shifted treatment values `A_t^d`.
    apply_policy_t = function(t, A_vec, H_df) {
      self$policies[[t]]$apply_policy(A_vec, H_df)
    },
    #' @description Compute the post-policy density `g_t^d(A_t | H_t)` from
    #'   a supplied conditional density estimator.
    #' @param t Time index (1-based).
    #' @param A_vec Vector of treatment values.
    #' @param H_df Data frame of history covariates.
    #' @param density_fun Function `(A_vec, H_df) -> g_t(A_t | H_t)`.
    #' @return Numeric vector of `g_t^d(A_t | H_t)`.
    gd_from_density_t = function(t, A_vec, H_df, density_fun) {
      self$policies[[t]]$gd_from_density(A_vec, H_df, density_fun = density_fun)
    }
  )
)

#' Repeat a single-time MTP across all time points
#'
#' @param mtp A single-time MTP object implementing `apply_policy` and
#'   `gd_from_density`.
#' @param tau Integer number of time points.
#'
#' @return An `LMTPPolicySequence` where the same MTP is applied at each time.
#' @export
repeat_policy_over_time <- function(mtp, tau) {
  LMTPPolicySequence$new(rep(list(mtp), tau))
}

#' Construct per-time MTPs from a factory
#'
#' @description
#' Build a list of per-time MTP objects by repeatedly calling a user-supplied
#' factory function that maps `(t, A_name, H_names)` to a new MTP instance.
#'
#' @param A_cols Character vector of treatment column names `c("A1", ..., "AT")`.
#' @param L_cols List of length `T`, each element a character vector of
#'   time-varying covariate names at that time.
#' @param W_cols Optional character vector of baseline covariate names.
#' @param mtp_factory Function with signature `(t, A_name, H_names)` returning
#'   a new MTP instance.
#'
#' @return An `LMTPPolicySequence`.
#' @export
make_per_time_policies <- function(
    A_cols,
    L_cols,
    W_cols = NULL,
    mtp_factory
) {
  stopifnot(length(A_cols) == length(L_cols))
  tau <- length(A_cols)
  policies <- vector("list", tau)
  for (t in seq_len(tau)) {
    A_name  <- A_cols[[t]]
    H_names <- c(W_cols %||% character(0),
                 unlist(L_cols[seq_len(t)], use.names = FALSE),
                 if (t > 1) A_cols[seq_len(t - 1)] else character(0))
    policies[[t]] <- mtp_factory(t, A_name, H_names)
  }
  LMTPPolicySequence$new(policies)
}


# LMTP Data Structure ----------------------------------------------------------

#' Data structure for longitudinal LMTP analyses
#'
#' @description
#' `LMTP_Data_Struct` is a light-weight R6 container for longitudinal data
#' organized as baseline covariates `W`, time-varying covariates `L_t`,
#' treatments `A_t`, and a final outcome `Y`. It provides convenient methods
#' to extract per-time histories `H_t`, treatments, and outcomes needed by
#' the LMTP estimators.
#'
#' @export
LMTP_Data_Struct <- R6::R6Class(
  "LMTP_Data_Struct",
  public = list(
    #' @field data Data frame containing `id`, `W`, `L_t`, `A_t`, and `Y`.
    data = NULL,
    #' @field id_col Name of the id column.
    id_col = NULL,
    #' @field n_timesteps Number of treatment times (tau).
    n_timesteps = NULL,
    #' @field A_cols Character vector of treatment column names.
    A_cols = NULL,
    #' @field L_cols List of length `tau` of time-varying covariate names.
    L_cols = NULL,
    #' @field W_cols Optional character vector of baseline covariate names.
    W_cols = NULL,
    #' @field Y_col Outcome column name.
    Y_col = NULL,

    #' @description Construct a new LMTP data structure.
    #' @param data Data frame with all variables.
    #' @param id_col Name of id column.
    #' @param n_timesteps Integer number of treatment times.
    #' @param A_cols Character vector of treatment column names.
    #' @param L_cols List of character vectors for time-varying covariates.
    #' @param W_cols Optional character vector of baseline covariate names.
    #' @param Y_col Outcome column name.
    initialize = function(data, id_col, n_timesteps, A_cols, L_cols,
                          W_cols = NULL, Y_col) {
      self$data <- data
      self$id_col <- id_col
      self$n_timesteps <- as.integer(n_timesteps)
      self$A_cols <- A_cols
      self$L_cols <- L_cols
      self$W_cols <- W_cols
      self$Y_col <- Y_col
      private$validate()
    },

    #' @description Number of time points.
    #' @return Integer `tau`.
    tau = function() self$n_timesteps,

    #' @description Get the history `H_t = (W, L_1..L_t, A_1..A_{t-1})`.
    #' @param t Time index (1-based).
    #' @return Data frame of covariates in the history.
    H = function(t) {
      stopifnot(t >= 1, t <= self$tau())
      A_hist <- if (t > 1) self$A_cols[seq_len(t - 1)] else character(0)
      L_hist <- unlist(self$L_cols[seq_len(t)], use.names = FALSE)
      cols <- c(A_hist, L_hist, if (!is.null(self$W_cols)) self$W_cols)
      self$data[, cols, drop = FALSE]
    },

    #' @description Get `(A_t, H_t)` as a single data frame.
    #' @param t Time index (1-based).
    #' @return Data frame with `H_t` and the treatment column `A_t`.
    AH = function(t) {
      df <- self$H(t)
      df[[self$A_cols[[t]]]] <- self$data[[self$A_cols[[t]]]]
      df
    },

    #' @description Get the treatment vector at time `t`.
    #' @param t Time index.
    #' @return Vector `A_t`.
    A = function(t) self$data[[self$A_cols[[t]]]],

    #' @description Get the time-varying covariates at time `t`.
    #' @param t Time index.
    #' @return Data frame of covariates `L_t`.
    L = function(t) self$data[, self$L_cols[[t]], drop = FALSE],

    #' @description Get baseline covariates `W`.
    #' @return Data frame of baseline covariates or `NULL`.
    W = function() if (is.null(self$W_cols)) NULL else self$data[, self$W_cols, drop = FALSE],

    #' @description Get the outcome vector `Y`.
    #' @return Outcome vector.
    Y = function() self$data[[self$Y_col]],

    #' @description Convenience: return a list with `A`, `H`, and `H_next` at time `t`.
    #' @param t Time index.
    #' @return List with elements `A`, `H`, `AH`, and `H_next`.
    view_t = function(t) {
      list(
        A = self$A(t),
        H = self$H(t),
        AH = self$AH(t),
        H_next = if (t < self$tau()) self$H(t + 1) else NULL
      )
    }
  ),
  private = list(
    validate = function() {
      stopifnot(self$id_col %in% colnames(self$data))
      stopifnot(length(self$A_cols) == self$n_timesteps)
      stopifnot(is.list(self$L_cols), length(self$L_cols) == self$n_timesteps)
      stopifnot(all(vapply(self$L_cols, function(v) is.character(v) && length(v) >= 1, TRUE)))
      needed <- c(unlist(self$L_cols), self$A_cols,
                  if (!is.null(self$W_cols)) self$W_cols,
                  self$Y_col)
      if (!all(needed %in% colnames(self$data))) {
        miss <- setdiff(needed, colnames(self$data))
        stop("LMTP_Data_Struct: missing columns in data: ", paste(miss, collapse = ", "))
      }
    }
  )
)

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

# TMLE fluctuation models ------------------------------------------------------

#' Identity-link TMLE fluctuation for LMTP
#'
#' @description
#' Implements a simple linear fluctuation on the identity scale. This is
#' mainly used for continuous outcomes; for bounded outcomes, consider
#' `LMTPFluctuationLogit`.
#'
#' @export
LMTPFluctuationIdentity <- R6::R6Class(
  "LMTPFluctuationIdentity",
  public = list(
    #' @description Fit a single-step fluctuation at time `t`.
    #' @param Q0_fun_t Function `(A,H) -> Q0_t(A,H)`.
    #' @param Qvec Numeric vector of `Q0_t(A_t, H_t)` at the observed data.
    #' @param target_vec Numeric target vector (e.g., `Y` or next-step
    #'   pseudo-outcome).
    #' @param K_obs_t Numeric clever covariate vector at the observed data.
    #' @param Kprov K-provider (unused here but kept for API compatibility).
    #' @param t Time index.
    #' @return A list with `epsilon` and `wrap(A,H)` returning `Q_t^*(A,H)`.
    fit_update = function(Q0_fun_t, Qvec, target_vec, K_obs_t, Kprov, t) {
      df <- data.frame(target = target_vec, offset = Qvec, K = K_obs_t)
      fit <- stats::glm(target ~ -1 + offset(offset) + K, data = df, family = gaussian())
      eps <- tryCatch(unname(coef(fit)["K"]), error = function(e) NA_real_)
      if (is.na(eps)) eps <- 0

      K_vec <- K_obs_t
      list(
        epsilon = eps,
        wrap = function(A_vec, H_df) {
          Q0 <- Q0_fun_t(A_vec, H_df)
          Q0 + eps * K_vec
        }
      )
    }
  )
)

#' Logit-link TMLE fluctuation for bounded LMTP outcomes
#'
#' @description
#' Fluctuation model on the logit scale with optional truncation and
#' bounds, appropriate for `[0, 1]`-bounded outcomes.
#'
#' @export
LMTPFluctuationLogit <- R6::R6Class(
  "LMTPFluctuationLogit",
  public = list(
    a = 0, b = 1, trunc = 1e-6,
    #' @description Create a new logit fluctuation object.
    #' @param bounds Two-element numeric vector giving the lower and upper
    #'   bounds of the outcome.
    #' @param truncation_alpha Small positive number used to truncate
    #'   probabilities away from 0 and 1.
    initialize = function(bounds = c(0, 1), truncation_alpha = 1e-6) {
      self$a <- bounds[1]; self$b <- bounds[2]; self$trunc <- truncation_alpha
    },
    .to01 = function(x) (x - self$a) / (self$b - self$a),
    .from01 = function(z) self$a + (self$b - self$a) * z,
    .b01 = function(z) pmin(pmax(z, self$trunc), 1 - self$trunc),

    #' @description Fit a single-step fluctuation at time `t`.
    #' @inherit LMTPFluctuationIdentity$fit_update params return
    fit_update = function(Q0_fun_t, Qvec, target_vec, K_obs_t, Kprov, t) {
      Y01 <- self$.to01(target_vec)
      off <- qlogis(self$.b01(self$.to01(Qvec)))
      df <- data.frame(Y01 = Y01, offset = off, K = K_obs_t)
      fit <- stats::glm(Y01 ~ -1 + offset(offset) + K, data = df, family = binomial())
      eps <- tryCatch(unname(coef(fit)["K"]), error = function(e) NA_real_)
      if (is.na(eps)) eps <- 0
      K_vec <- K_obs_t
      list(
        epsilon = eps,
        wrap = function(A_vec, H_df) {
          Q0 <- Q0_fun_t(A_vec, H_df)
          z  <- qlogis(self$.b01(self$.to01(Q0))) + eps * K_vec
          self$.from01(plogis(z))
        }
      )
    }
  )
)

# TMLE / G-computation / IPW for LMTP -----------------------------------------

#' LMTP estimation via TMLE, outcome regression, or IPW
#'
#' @description
#' Core workhorse to estimate the marginal LMTP parameter
#' `psi = E[Y^d]` under a longitudinal modified treatment policy. This
#' function performs a sequential outcome regression, optionally applies
#' a TMLE fluctuation using the LMTP clever covariate, and computes the
#' efficient influence curve using the *unfluctuated* outcome regression
#' trajectory `Q0_t` (as recommended for variance estimation in
#' longitudinal TMLE; see discussion in Díaz et al. 2023 and related
#' work). The plug-in estimate is still based on the fluctuated
#' trajectory when `method = "tmle"`.
#'
#' @param ds `LMTP_Data_Struct` instance.
#' @param policy_seq `LMTPPolicySequence` describing the LMTP.
#' @param learners_Q List (or single spec) of learner specifications for
#'   sequential outcome regressions.
#' @param learners_g_factory An `LMTPNuisanceFactory` instance (already
#'   configured, but not yet trained).
#' @param fml_Q Optional formula or list of formulas for the outcome
#'   regressions.
#' @param outcome_link Either `"identity"` or `"logit"`.
#' @param bounds Numeric vector of length 2 giving outcome bounds when
#'   `outcome_link = "logit"`.
#' @param maxit Maximum number of outer TMLE iterations.
#' @param eps_tol Convergence tolerance for the fluctuation (ignored when
#'   `maxit = 1`).
#' @param repeat_lnrs Logical; if `TRUE`, recycle `learners_Q` across
#'   time points when only one is supplied.
#' @param method One of `"tmle"`, `"outcome"` (sequential G-computation
#'   without targeting), or `"ipw"`.
#'
#' @return An object of class `tmle_lmtp` with components:
#'   * `psi`: point estimate of `E[Y^d]`.
#'   * `se`: standard error.
#'   * `ci95`: 95\% Wald confidence interval.
#'   * `ic`: estimated influence curve values.
#'   * `Q_star`: list of fluctuated regression functions `Q_t^*`.
#'   * `Q0`: list of unfluctuated regression functions `Q0_t`.
#'   * `Kprov`: the K-provider used.
#'
#' @export
fit_tmle_for_LMTP <- function(
    ds,
    policy_seq,
    learners_Q,
    learners_g_factory,
    fml_Q = NULL,
    outcome_link = c("identity", "logit"),
    bounds = c(0, 1),
    maxit = 1,
    eps_tol = 1e-6,
    repeat_lnrs = TRUE,
    method = c("tmle", "outcome", "ipw")
) {
  method <- match.arg(method)
  outcome_link <- match.arg(outcome_link)
  tau <- ds$tau()

  # 1) Fit g and r once
  learners_g_factory$train(ds)

  # 2) Helpers from the factory
  Kprov <- learners_g_factory$k_provider(ds)
  learners_Q <- .rep_if_needed(learners_Q, ds$tau(), repeat_lnrs)

  Q_trainer_fun <- learners_g_factory$q_trainer(
    learners_Q = if (length(learners_Q) == 1) rep(list(learners_Q[[1]]), tau) else learners_Q,
    fml_Q = fml_Q,
    repeat_fmls = is.null(fml_Q) || length(fml_Q) == 1
  )

  fluct <- switch(outcome_link,
                  identity = LMTPFluctuationIdentity$new(),
                  logit    = LMTPFluctuationLogit$new(bounds = bounds))

  # 3) Initialize storage for Q_star and Q0 closures
  Q_star <- vector("list", tau)
  Q0_list <- vector("list", tau)

  # 4) Outer iterate
  iter <- 1; eps <- 1e6
  while (iter <= maxit && eps > eps_tol) {

    tilde_next <- ds$Y()

    for (t in tau:1) {
      Ht <- ds$H(t)

      # Train initial Q0_t on pseudo-outcome
      Q0_t <- Q_trainer_fun(ds, t, tilde_next)
      Q0_list[[t]] <- Q0_t

      # target for epsilon fit
      target_vec <- if (t < tau) {
        Hnext <- ds$H(t + 1)
        Astar_next <- policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext)
        if (!is.null(Q_star[[t + 1]])) {
          Q_star[[t + 1]](Astar_next, Hnext)
        } else {
          Q0_list[[t + 1]](Astar_next, Hnext)
        }
      } else {
        ds$Y()
      }

      # clever covariate and offset evals on observed data
      K_obs_t <- Kprov$K_obs(t)
      Qvec <- Q0_t(ds$A(t), Ht)

      if (method == "tmle") {
        up <- fluct$fit_update(Q0_fun_t = Q0_t,
                               Qvec = Qvec,
                               target_vec = target_vec,
                               K_obs_t = K_obs_t,
                               Kprov = Kprov, t = t)
        eps <- up$epsilon
        Q_star[[t]] <- up$wrap
      } else {
        Q_star[[t]] <- Q0_t
        eps <- 0
      }

      # carry upstream the pseudo-outcome under the policy
      Astar_t <- policy_seq$apply_policy_t(t, ds$A(t), Ht)
      tilde_next <- Q_star[[t]](Astar_t, Ht)
    }

    iter <- iter + 1
  }

  # 5) Estimate psi
  H1 <- ds$H(1)
  Astar1 <- policy_seq$apply_policy_t(1, ds$A(1), H1)
  if (method %in% c("tmle", "outcome")) {
    psi_hat <- mean(Q_star[[1]](Astar1, H1))
  } else if (method == "ipw") {
    psi_hat <- mean(Kprov$K_obs(tau) * ds$Y())
  }

  # 6) Influence curve using *unfluctuated* Q0_t
  if (method == "tmle") {
    n <- nrow(ds$data)
    phi <- rep(0, n)
    for (t in seq_len(tau)) {
      Ht <- ds$H(t)
      mt_obs <- Q0_list[[t]](ds$A(t), Ht)
      mnext <- if (t < tau) {
        Hnext <- ds$H(t + 1)
        Astar_next <- policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext)
        Q0_list[[t + 1]](Astar_next, Hnext)
      } else ds$Y()
      w <- Kprov$K_obs(t)
      phi <- phi + w * (mnext - mt_obs)
    }
    plug <- Q_star[[1]](Astar1, H1)
    ic <- phi + plug - psi_hat
    se <- sqrt(stats::var(ic) / n)
    ci95 <- c(psi_hat - 1.96 * se, psi_hat + 1.96 * se)
  } else {
    se <- NA; ci95 <- c(NA, NA); ic <- NA
  }

  structure(list(
    psi = psi_hat,
    se = se,
    ci95 = ci95,
    ic = ic,
    Q_star = Q_star,
    Q0 = Q0_list,
    Kprov = Kprov
  ), class = c("tmle_lmtp", "list"))
}

# Heterogeneous LMTP effects ---------------------------------------------------

#' Heterogeneous LMTP effects for baseline subgroups
#'
#' @description
#' Estimate subgroup-specific LMTP effects
#'  \eqn{\theta_v = E[Y^d | V = v]} for a set of baseline subgroups defined by
#' baseline covariates `V ⊂ W`. The procedure implements the joint TMLE
#' described in the prompt: a single set of nuisance fits and a
#' multi-dimensional clever covariate that targets all subgroup effects
#' simultaneously, following the strategy of Wei et al. (2023) adapted to
#' the LMTP setting.
#'
#' Subgroups are defined as the unique combinations of columns in `V`.
#' For each subgroup level `v`, the target parameter is
#' `theta_v = E[Y^d | V = v]` and the EIF is
#'
#' \deqn{D_v^*(O) = \frac{I(V=v)}{p_v}\left\{\sum_{t=1}^\tau
#'   \omega_t( m_{t+1}(A_t^d, H_t) - E[m_{t+1}(A_{t+1}, H_{t+1}) | A_t, H_t]) + m_1(A_1^d, L_1) - \theta_v \right\},}
#'
#' with \eqn{p_v = P(V = v)} and \eqn{\omega_t = \prod_{s=1}^t r_s(A_s, H_s)}.
#' In this implementation, the TMLE produces a targeted (fluctuated) substitution estimator \eqn{m_1^*},
#' and uses an unfluctuated substitution estimators \eqn{m_1} in the influence curve used for variance
#' calculations and inference.
#'
#' @param ds `LMTP_Data_Struct` instance.
#' @param policy_seq `LMTPPolicySequence` describing the LMTP.
#' @param learners_Q List (or single spec) of learner specifications for
#'   the outcome regressions.
#' @param learners_g_factory An `LMTPNuisanceFactory` instance (already
#'   configured, but not yet trained).
#' @param V Character vector of baseline column names in `ds$W()` used to
#'   define subgroups. All unique combinations of these columns define
#'   the set \eqn{\mathcal V} over which subgroup effects are estimated.
#' @param fml_Q Optional formula or list of formulas for the outcome
#'   regressions.
#' @param outcome_link Either `"identity"` or `"logit"`.
#' @param bounds Bounds for the outcome when `outcome_link = "logit"`.
#' @param maxit Currently ignored (one-step TMLE); included for future
#'   compatibility.
#' @param eps_tol Currently ignored; included for compatibility.
#' @param repeat_lnrs Logical; if `TRUE`, recycle `learners_Q` across
#'   time points when only one is supplied.
#' @param method One of `"tmle"`, `"outcome"`, or `"ipw"`. When
#'   `method = "ipw"`, only inverse probability weighted estimators are
#'   computed, without outcome regression.
#' @param stabilization Character string controlling potential future
#'   stabilization of the subgroup clever covariates / weights. Currently
#'   only `"naive"` is implemented; `"dual"` is reserved for a future
#'   implementation of the Lagrangian-dual approach of Wei et al. (2023).
#'
#' @return An object of class `tmle_lmtp_heterogeneous` with components:
#'   * `theta`: vector of subgroup-specific estimates.
#'   * `se`: standard errors for each subgroup.
#'   * `ci95`: matrix of 95\% Wald intervals (columns `lower`, `upper`).
#'   * `ic`: `n x J` matrix of influence curve values, with columns
#'      corresponding to subgroups.
#'   * `V_levels`: factor levels of the subgroup variable (unique
#'      combinations of `V`).
#'   * `p_v`: empirical subgroup proportions.
#'   * `method`: estimation method used.
#'   * `subgroup_vars`: names of the baseline variables that define the
#'     subgroups.
#'
#' @export
fit_lmtp_heterogeneous <- function(
    ds,
    policy_seq,
    learners_Q,
    learners_g_factory,
    V,
    fml_Q = NULL,
    outcome_link = c("identity", "logit"),
    bounds = c(0, 1),
    maxit = 1,
    eps_tol = 1e-6,
    repeat_lnrs = TRUE,
    method = c("tmle", "outcome", "ipw"),
    stabilization = c("naive", "dual")) {

  outcome_link <- match.arg(outcome_link)
  method <- match.arg(method)
  stabilization <- match.arg(stabilization)

  if (stabilization != "naive") {
    stop("stabilization = 'dual' is not yet implemented; only 'naive' is currently supported.")
  }

  tau <- ds$tau()
  n <- nrow(ds$data)

  # Subgroup variable V: factor of unique baseline combinations
  W_df <- ds$W()
  if (is.null(W_df)) stop("No baseline W_cols supplied in LMTP_Data_Struct; cannot define V.")
  if (!all(V %in% colnames(W_df))) {
    missing_V <- setdiff(V, colnames(W_df))
    stop("Columns not found in baseline W: ", paste(missing_V, collapse = ", "))
  }
  V_df <- W_df[, V, drop = FALSE]
  V_str <- interaction(V_df, drop = TRUE, lex.order = TRUE)
  V_levels <- levels(V_str)
  J <- length(V_levels)

  p_v <- prop.table(table(V_str))
  p_v_vec <- as.numeric(p_v)[match(V_levels, names(p_v))]

  # Train nuisance functions
  learners_g_factory$train(ds)
  Kprov <- learners_g_factory$k_provider(ds)
  learners_Q <- .rep_if_needed(learners_Q, ds$tau(), repeat_lnrs)

  Q_trainer_fun <- learners_g_factory$q_trainer(
    learners_Q = if (length(learners_Q) == 1) rep(list(learners_Q[[1]]), tau) else learners_Q,
    fml_Q = fml_Q,
    repeat_fmls = is.null(fml_Q) || length(fml_Q) == 1
  )

  # Precompute omega_t and subgroup clever covariates K_{t,v}
  omega_list <- lapply(seq_len(tau), function(t) Kprov$K_obs(t))
  K_list <- vector("list", tau)
  for (t in seq_len(tau)) {
    omega_t <- omega_list[[t]]
    K_t <- matrix(0, nrow = n, ncol = J)
    colnames(K_t) <- paste0("K", seq_len(J))
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      if (length(idx) > 0) {
        K_t[idx, j] <- omega_t[idx] / p_v_vec[j]
      }
    }
    K_list[[t]] <- K_t
  }

  # Storage for Q0 and Q* trajectories (vectors evaluated on the sample)
  Q0_obs_list <- vector("list", tau)
  Q0_d_list   <- vector("list", tau)
  Qstar_obs_list <- vector("list", tau)
  Qstar_d_list   <- vector("list", tau)

  # One-step TMLE / G-computation backward recursion
  tilde_next <- ds$Y()

  for (t in tau:1) {
    Ht <- ds$H(t)
    At <- ds$A(t)
    Astar_t <- policy_seq$apply_policy_t(t, At, Ht)

    Q0_t <- Q_trainer_fun(ds, t, tilde_next)
    Q0_obs <- Q0_t(At, Ht)
    Q0_d   <- Q0_t(Astar_t, Ht)

    Q0_obs_list[[t]] <- Q0_obs
    Q0_d_list[[t]]   <- Q0_d

    if (method == "ipw") {
      # no targeting / outcome regression used only for diagnostics
      Qstar_obs <- Q0_obs
      Qstar_d   <- Q0_d
    } else {
      K_t <- K_list[[t]]
      target_vec <- tilde_next

      if (outcome_link == "identity") {
        df <- data.frame(target = target_vec,
                         offset = Q0_obs,
                         K_t)
        fit <- stats::glm(target ~ -1 + offset(offset) + ., data = df, family = gaussian())
        eps <- coef(fit)
        eps[is.na(eps)] <- 0
        eps <- eps[setdiff(names(eps), "offset")]  # keep only K coefficients
        eps_vec <- as.numeric(eps)
        eps_vec <- if (length(eps_vec) < J) c(eps_vec, rep(0, J - length(eps_vec))) else eps_vec
        Qstar_obs <- as.numeric(Q0_obs + K_t %*% eps_vec)
        Qstar_d   <- as.numeric(Q0_d   + K_t %*% eps_vec)
      } else {
        # logit link with bounds
        a <- bounds[1]; b <- bounds[2]
        to01 <- function(x) (x - a) / (b - a)
        from01 <- function(z) a + (b - a) * z
        b01 <- function(z) pmin(pmax(z, 1e-6), 1 - 1e-6)

        Y01 <- to01(target_vec)
        off <- qlogis(b01(to01(Q0_obs)))
        df <- data.frame(Y01 = Y01,
                         offset = off,
                         K_t)
        fit <- stats::glm(Y01 ~ -1 + offset(offset) + ., data = df, family = binomial())
        eps <- coef(fit)
        eps[is.na(eps)] <- 0
        eps <- eps[setdiff(names(eps), "offset")]  # only K coefficients
        eps_vec <- as.numeric(eps)
        eps_vec <- if (length(eps_vec) < J) c(eps_vec, rep(0, J - length(eps_vec))) else eps_vec

        z_obs <- off + as.numeric(K_t %*% eps_vec)
        z_d   <- qlogis(b01(to01(Q0_d))) + as.numeric(K_t %*% eps_vec)
        Qstar_obs <- from01(plogis(z_obs))
        Qstar_d   <- from01(plogis(z_d))
      }
    }

    Qstar_obs_list[[t]] <- Qstar_obs
    Qstar_d_list[[t]]   <- Qstar_d

    # Update pseudo-outcome for next step
    tilde_next <- if (method == "tmle") Qstar_d else Q0_d
  }

  # Plug-in estimators for each subgroup
  theta <- numeric(J)
  names(theta) <- V_levels

  if (method %in% c("tmle", "outcome")) {
    Q1_d <- if (method == "tmle") Qstar_d_list[[1]] else Q0_d_list[[1]]
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      theta[j] <- mean(Q1_d[idx] / p_v_vec[j])
    }
  } else if (method == "ipw") {
    omega_tau <- omega_list[[tau]]
    Y <- ds$Y()
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      theta[j] <- mean(omega_tau[idx] * Y[idx] / p_v_vec[j])
    }
  }

  # Influence curve for tmle / outcome using unfluctuated Q0_t
  if (method %in% c("tmle", "outcome")) {
    phi <- rep(0, n)
    for (t in seq_len(tau)) {
      mt_obs <- Q0_obs_list[[t]]
      mnext <- if (t < tau) Q0_d_list[[t + 1]] else ds$Y()
      phi <- phi + omega_list[[t]] * (mnext - mt_obs)
    }
    plug <- Q0_d_list[[1]]
    ic_mat <- matrix(0, nrow = n, ncol = J)
    colnames(ic_mat) <- V_levels
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      ic_v <- phi + plug - theta[j]
      ic_mat[, j] <- 0
      ic_mat[idx, j] <- ic_v[idx] / p_v_vec[j]
    }
    se <- sqrt(colSums((ic_mat - matrix(colMeans(ic_mat), nrow = n, ncol = J, byrow = TRUE))^2) / (n - 1) / n)
    ci95 <- cbind(lower = theta - 1.96 * se,
                  upper = theta + 1.96 * se)
  } else {
    ic_mat <- matrix(NA_real_, nrow = n, ncol = J)
    colnames(ic_mat) <- V_levels
    se <- rep(NA_real_, J)
    ci95 <- cbind(lower = rep(NA_real_, J), upper = rep(NA_real_, J))
  }

  res <- list(
    theta = theta,
    se = se,
    ci95 = ci95,
    ic = ic_mat,
    V_levels = V_levels,
    p_v = p_v_vec,
    method = method,
    subgroup_vars = V
  )
  class(res) <- c("tmle_lmtp_heterogeneous", "list")
  res
}
