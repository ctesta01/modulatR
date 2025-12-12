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
