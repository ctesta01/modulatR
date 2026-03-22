#' Scalar parameter object for LMTP estimation
#'
#' @description
#' `LMTPParameter` is an R6 class representing a scalar target parameter in a
#' longitudinal modified treatment policy problem.
#'
#' The class is intentionally narrow. It does not fit nuisance functions or
#' perform targeting itself. Instead, it defines:
#' \itemize{
#'   \item the parameter name,
#'   \item the scale of the target parameter,
#'   \item how to compute the plug-in estimate from fitted regressions,
#'   \item how to compute the efficient influence function contribution once
#'   nuisance functions are available.
#' }
#'
#' The default implementation corresponds to the LMTP mean parameter
#' \eqn{\psi = \mathbb{E}[Y^d]}.
#'
#' @export
#' @examples
#' set.seed(1)
#' n <- 200
#'
#' df <- data.frame(
#'   W1 = rnorm(n),
#'   L1 = rnorm(n)
#' )
#' df$A1 <- rbinom(n, 1, plogis(0.5 * df$W1 - 0.7 * df$L1))
#' df$L2 <- rnorm(n, mean = 0.3 * df$A1 + 0.2 * df$W1)
#' df$A2 <- rbinom(n, 1, plogis(-0.3 * df$W1 + 0.6 * df$L2))
#' df$Y  <- rnorm(n, mean = 0.5 * df$A1 + 0.7 * df$A2 + 0.3 * df$L2)
#'
#' ds <- LMTPData$new(
#'   data = df,
#'   A_cols = c("A1", "A2"),
#'   L_cols = list("L1", "L2"),
#'   W_cols = "W1",
#'   Y_col = "Y"
#' )
#'
#' pol <- mtp_discrete(
#'   map_fun = function(A, H) 1 - A,
#'   support = c(0, 1)
#' )
#' policy_seq <- repeat_policy_over_time(pol, tau = ds$tau())
#'
#' nf <- LMTPNuisanceFactory$new(
#'   learners_g = nadir::lnr_logistic,
#'   policy_seq = policy_seq,
#'   A_type = "discrete",
#'   g_mode = "density"
#' )
#' nf$train(ds)
#'
#' q_fit_factory <- nf$q_trainer(
#'   learners_Q = nadir::lnr_glm
#' )
#'
#' Q2 <- q_fit_factory(ds, t = 2, pseudo_outcome_vec = ds$Y())
#' Q1 <- q_fit_factory(ds, t = 1, pseudo_outcome_vec = Q2(ds$A(2), ds$H(2)))
#'
#' param <- LMTPMean(outcome_type = "continuous")
#' psi <- param$plugin(ds, policy_seq, Q_list = list(Q1, Q2))
#' psi
#'
#' kp <- nf$k_provider(ds)
#' head(param$eif(ds, policy_seq, Q_list = list(Q1, Q2), k_provider = kp))
LMTPParameter <- R6::R6Class(
  "LMTPParameter",
  public = list(
    #' @field name Character name of the parameter.
    name = NULL,

    #' @field outcome_type Character outcome type, typically `"continuous"` or `"binary"`.
    outcome_type = NULL,

    #' @field scale Character scale on which the parameter lives.
    #'   Usually `"identity"`.
    scale = NULL,

    #' @description Create a new LMTP parameter object.
    #' @param name Parameter name.
    #' @param outcome_type Outcome type.
    #' @param scale Parameter scale.
    initialize = function(name = "E[Y^d]",
                          outcome_type = c("continuous", "binary"),
                          scale = c("identity")) {
      self$name <- name
      self$outcome_type <- match.arg(outcome_type)
      self$scale <- match.arg(scale)
      invisible(self)
    },

    #' @description Number of treatment times.
    #' @param ds An `LMTPData` object.
    #' @return Integer `tau`.
    tau = function(ds) {
      if (!inherits(ds, "LMTPData")) {
        stop("`ds` must inherit from `LMTPData`.")
      }
      ds$tau()
    },

    #' @description Compute the plug-in estimate from the final-stage regression.
    #'
    #' @details
    #' For the LMTP mean parameter, the plug-in estimate is the empirical mean
    #' of the stage-1 shifted regression:
    #' \deqn{
    #'   \psi_n = \frac{1}{n}\sum_{i=1}^n \bar Q_1^d(H_{1,i}).
    #' }
    #'
    #' This method assumes `Q_star_list[[1]]` is a function of `(A, H)` returning
    #' the shifted regression at time 1 evaluated at `A = d_1(A_1, H_1)`.
    #'
    #' @param ds An `LMTPData` object.
    #' @param policy_seq An `LMTPPolicySequence`.
    #' @param Q_list List of fitted outcome regression functions. Each
    #'   `Q_list[[t]]` should be a function `(A_vec, H_df) -> numeric`.
    #' @return Numeric scalar estimate.
    plugin = function(ds, policy_seq, Q_list) {
      private$validate_common_inputs(ds, policy_seq, Q_list)

      A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))
      Q1_star <- Q_list[[1]](A1_star, ds$H(1))
      mean(Q1_star)
    },

    #' @description Compute the un-targeted EIF contribution.
    #'
    #' @details
    #' For the LMTP mean parameter, the EIF has the recursive form
    #' \deqn{
    #' D(O) =
    #' \sum_{t=1}^{\tau} \omega_t \left\{\bar Q_{t+1}^d - \bar Q_t^d\right\}
    #' + \bar Q_1^d - \psi,
    #' }
    #' where by convention \eqn{\bar Q_{\tau+1}^d = Y}.
    #'
    #' This method expects:
    #' \itemize{
    #'   \item `Q_list[[t]]` gives \eqn{\bar Q_t^d(A_t, H_t)} via a function of `(A, H)`,
    #'   \item `k_provider$K_obs(t)` returns \eqn{\omega_t},
    #'   \item `policy_seq$apply_t(t, ...)` computes the shifted treatment at time `t`.
    #' }
    #'
    #' @param ds An `LMTPData` object.
    #' @param policy_seq An `LMTPPolicySequence`.
    #' @param Q_list List of fitted outcome regression functions.
    #' @param k_provider Output of `LMTPNuisanceFactory$k_provider(ds)`.
    #' @param psi Optional scalar plug-in estimate. If `NULL`, it is computed
    #'   internally via `$plugin()`.
    #' @return Numeric vector of length `n`.
    eif = function(ds, policy_seq, Q_list, k_provider, psi = NULL) {
      private$validate_common_inputs(ds, policy_seq, Q_list)

      if (is.null(k_provider) || !is.function(k_provider$K_obs)) {
        stop("`k_provider` must contain a function `K_obs(t)`.")
      }

      if (is.null(psi)) {
        psi <- self$plugin(ds, policy_seq, Q_list)
      }

      tau <- ds$tau()
      D <- rep(0, ds$n)

      for (t in seq_len(tau)) {
        A_star_t <- policy_seq$apply_t(t, ds$A(t), ds$H(t))
        Q_t_star <- Q_list[[t]](A_star_t, ds$H(t))

        if (t < tau) {
          A_star_next <- policy_seq$apply_t(t + 1, ds$A(t + 1), ds$H(t + 1))
          Q_next_star <- Q_list[[t + 1]](A_star_next, ds$H(t + 1))
        } else {
          Q_next_star <- ds$Y()
        }

        D <- D + k_provider$K_obs(t) * (Q_next_star - Q_t_star)
      }

      A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))
      Q1_star <- Q_list[[1]](A1_star, ds$H(1))

      D + Q1_star - psi
    },

    #' @description Compute estimate, EIF, standard error, and Wald CI.
    #' @param ds An `LMTPData` object.
    #' @param policy_seq An `LMTPPolicySequence`.
    #' @param Q_list List of fitted outcome regression functions.
    #' @param k_provider Output of `LMTPNuisanceFactory$k_provider(ds)`.
    #' @param alpha Significance level for confidence interval.
    #' @return A named list with `estimate`, `eif`, `std_error`, `var`, and `ci`.
    infer = function(ds,
                     policy_seq,
                     Q_list,
                     k_provider,
                     alpha = 0.05) {
      psi <- self$plugin(ds, policy_seq, Q_list)
      eif <- self$eif(ds, policy_seq, Q_list, k_provider, psi = psi)

      var_hat <- stats::var(eif) / ds$n
      se_hat <- sqrt(var_hat)
      z <- stats::qnorm(1 - alpha / 2)

      list(
        estimate = psi,
        eif = eif,
        var = var_hat,
        std_error = se_hat,
        ci = c(
          lower = psi - z * se_hat,
          upper = psi + z * se_hat
        )
      )
    },

    #' @description Print a compact summary.
    #' @return The object invisibly.
    print = function(...) {
      cat("LMTPParameter\n")
      cat("  name: ", self$name, "\n", sep = "")
      cat("  outcome_type: ", self$outcome_type, "\n", sep = "")
      cat("  scale: ", self$scale, "\n", sep = "")
      invisible(self)
    }
  ),

  private = list(
    validate_common_inputs = function(ds, policy_seq, Q_list) {
      if (!inherits(ds, "LMTPData")) {
        stop("`ds` must inherit from `LMTPData`.")
      }

      if (!inherits(policy_seq, "LMTPPolicySequence")) {
        stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
      }

      policy_seq$validate_against_data(ds)

      if (!is.list(Q_list) || length(Q_list) != ds$tau()) {
        stop("`Q_list` must be a list of length `tau`.")
      }

      ok <- vapply(Q_list, is.function, logical(1))
      if (!all(ok)) {
        stop("Each element of `Q_list` must be a function `(A_vec, H_df) -> numeric`.")
      }
    }
  )
)

#' LMTP mean parameter
#'
#' @description
#' Convenience constructor for the scalar LMTP mean parameter
#' \eqn{\psi = \mathbb{E}[Y^d]}.
#'
#' @param outcome_type Outcome type.
#' @return An `LMTPParameter` object.
#' @export
LMTPMean <- function(outcome_type = c("continuous", "binary")) {
  LMTPParameter$new(
    name = "E[Y^d]",
    outcome_type = match.arg(outcome_type),
    scale = "identity"
  )
}
