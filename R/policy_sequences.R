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
