
# Nuisance Factories for LMTP ------------------------------------------------

#' LMTP nuisance-function factory
#'
#' @description
#' A registry and adapter for nuisance functions used in
#' longitudinal modified treatment policy estimation.
#'
#' This object accepts black-box training functions. Each training function should
#' take one argument, an `LMTPData` object, and return predictions.
#'
#' The notation follows Díaz et al. (2023):
#'
#' - `m_t(a_t, h_t)` is the sequential outcome regression.
#' - `g_t(a_t | h_t)` is the observed treatment density or mass function.
#' - `g_t^d(a_t | h_t)` is the post-policy density or mass function.
#' - `r_t(a_t, h_t) = g_t^d(a_t | h_t) / g_t(a_t | h_t)` is the density ratio.
#' - `omega_t = prod_{s=1}^t r_s(A_s, H_s)` is the cumulative ratio.
#'
#' @export
LMTPNuisanceFactory <- R6::R6Class(
  "LMTPNuisanceFactory",
  public = list(
    policy_seq = NULL,

    # User-supplied learner functions
    r_fit_funs = NULL,
    m_fit_funs = NULL,

    # Numerical safeguards
    truncate_ratio = NULL,

    # Fitted ratio artifacts
    r_fits = NULL,
    r_obs = NULL,
    r_d = NULL,
    omega_obs = NULL,

    # Fitted outcome-regression artifacts
    m_fits = NULL,
    m_obs = NULL,
    m_d = NULL,

    initialize = function(policy_seq,
                          r_fit_funs,
                          m_fit_funs,
                          truncate_ratio = c(1e-6, Inf)) {
      self$policy_seq <- policy_seq
      self$r_fit_funs <- r_fit_funs
      self$m_fit_funs <- m_fit_funs
      self$truncate_ratio <- truncate_ratio

      private$validate()
      invisible(self)
    },

    train = function(data) {
      if (!inherits(data, "LMTPData")) {
        stop("`data` must inherit from `LMTPData`.")
      }

      self$policy_seq$validate_against_data(data)

      tau <- data$tau()
      r_fit_funs <- private$as_time_list(self$r_fit_funs, tau, "r_fit_funs")

      self$r_fits <- vector("list", tau)
      self$r_obs <- vector("list", tau)
      self$r_d <- vector("list", tau)
      self$omega_obs <- vector("list", tau)

      omega_prev <- rep(1, data$n)

      for (t in seq_len(tau)) {
        data_t <- data |>
          .with_policy_sequence(self$policy_seq) |>
          .with_lmtp_time(t)

        fit_t <- r_fit_funs[[t]](data_t)
        ratio_t <- private$extract_ratio_vectors(fit_t, t, data$n)

        r_obs_t <- .truncate_interval(
          ratio_t$r_obs,
          truncate = self$truncate_ratio
        )
        r_d_t <- .truncate_interval(
          ratio_t$r_d,
          truncate = self$truncate_ratio
        )

        omega_prev <- omega_prev * r_obs_t

        self$r_fits[[t]] <- fit_t
        self$r_obs[[t]] <- r_obs_t
        self$r_d[[t]] <- r_d_t
        self$omega_obs[[t]] <- omega_prev
      }

      invisible(self)
    },

    fit_m = function(data,
                     t,
                     pseudo_outcome,
                     pseudo_outcome_col = "..pseudo_outcome") {
      if (!inherits(data, "LMTPData")) {
        stop("`data` must inherit from `LMTPData`.")
      }
      if (is.null(self$m_fit_funs)) {
        stop("`m_fit_funs` is NULL.")
      }

      m_fit_funs <- private$as_time_list(
        self$m_fit_funs,
        data$tau(),
        "m_fit_funs"
      )

      data_t <- data |>
        .with_policy_sequence(self$policy_seq) |>
        .with_pseudo_outcome(
          pseudo_outcome = pseudo_outcome,
          pseudo_outcome_col = pseudo_outcome_col
        ) |>
        .with_lmtp_time(t)

      fit_t <- m_fit_funs[[t]](data_t)
      m_t <- private$extract_m_vectors(fit_t, t, data$n)

      if (is.null(self$m_fits)) self$m_fits <- vector("list", data$tau())
      if (is.null(self$m_obs)) self$m_obs <- vector("list", data$tau())
      if (is.null(self$m_d)) self$m_d <- vector("list", data$tau())

      self$m_fits[[t]] <- fit_t
      self$m_obs[[t]] <- m_t$m_obs
      self$m_d[[t]] <- m_t$m_d

      list(
        fit = fit_t,
        m_obs = m_t$m_obs,
        m_d = m_t$m_d
      )
    },

    omega = function(t) {
      private$check_trained_ratios()
      self$omega_obs[[t]]
    },

    omega_prev = function(t, n = NULL) {
      private$check_trained_ratios()

      if (t == 1L) {
        if (is.null(n)) {
          n <- length(self$omega_obs[[1L]])
        }
        return(rep(1, n))
      }

      self$omega_obs[[t - 1L]]
    },

    print = function(...) {
      cat("LMTPNuisanceFactory\n")
      cat("  stores: vector-valued r, omega, m\n", sep = "")
      cat("  ratios trained: ", !is.null(self$r_obs), "\n", sep = "")
      cat("  outcome regressions trained: ", !is.null(self$m_obs), "\n", sep = "")
      invisible(self)
    }
  ),

  private = list(
    validate = function() {
      if (!inherits(self$policy_seq, "LMTPPolicySequence")) {
        stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
      }

      if (is.null(self$r_fit_funs)) {
        stop("`r_fit_funs` must be supplied.")
      }

      if (is.null(self$m_fit_funs)) {
        stop("`m_fit_funs` must be supplied.")
      }

      if (!is.numeric(self$truncate_ratio) ||
          length(self$truncate_ratio) != 2L ||
          self$truncate_ratio[1] <= 0 ||
          self$truncate_ratio[2] < self$truncate_ratio[1]) {
        stop("`truncate_ratio` must be c(lower, upper), with 0 < lower <= upper.")
      }

      invisible(TRUE)
    },

    as_time_list = function(x, tau, label) {
      if (is.function(x)) {
        return(rep(list(x), tau))
      }

      if (is.list(x) &&
          length(x) == tau &&
          all(vapply(x, is.function, logical(1)))) {
        return(x)
      }

      stop("`", label, "` must be a function or a list of functions of length tau.")
    },

    extract_ratio_vectors = function(fit, t, n) {
      if (!is.list(fit)) {
        stop("Ratio learner at t = ", t, " must return a list.")
      }

      if (!is.null(fit$r_obs) && !is.null(fit$r_d)) {
        r_obs <- as.numeric(fit$r_obs)
        r_d <- as.numeric(fit$r_d)
      } else if (!is.null(fit$u_obs) && !is.null(fit$u_d)) {
        u_obs <- .clip_probability(as.numeric(fit$u_obs))
        u_d <- .clip_probability(as.numeric(fit$u_d))

        r_obs <- u_obs / (1 - u_obs)
        r_d <- u_d / (1 - u_d)
      } else {
        stop(
          "Ratio learner at t = ", t,
          " must return `r_obs` and `r_d`, or `u_obs` and `u_d`."
        )
      }

      if (length(r_obs) != n || length(r_d) != n) {
        stop(
          "Ratio learner at t = ", t,
          " must return vectors of length n."
        )
      }

      list(r_obs = r_obs, r_d = r_d)
    },

    extract_m_vectors = function(fit, t, n) {
      if (!is.list(fit)) {
        stop("Outcome learner at t = ", t, " must return a list.")
      }

      if (is.null(fit$m_obs) || is.null(fit$m_d)) {
        stop(
          "Outcome learner at t = ", t,
          " must return `m_obs` and `m_d`."
        )
      }

      m_obs <- as.numeric(fit$m_obs)
      m_d <- as.numeric(fit$m_d)

      if (length(m_obs) != n || length(m_d) != n) {
        stop(
          "Outcome learner at t = ", t,
          " must return vectors of length n."
        )
      }

      list(m_obs = m_obs, m_d = m_d)
    },

    check_trained_ratios = function() {
      if (is.null(self$r_obs) || is.null(self$r_d) || is.null(self$omega_obs)) {
        stop("Call `$train(data)` before requesting ratio or omega vectors.")
      }
      invisible(TRUE)
    }
  )
)



