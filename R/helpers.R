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

#' Use first argument if not null, else use second argument
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


.wald_ci <- function(est, se, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  c(lower = est - z * se, upper = est + z * se)
}

#' @export
identity_policy <- function(A_type = c("continuous", "discrete"),
                             discrete_support = NULL,
                             tau = 1L) {
  A_type <- match.arg(A_type)

  if (A_type == "continuous") {
    mtp_id <- MTP$new(
      treatment_type = "continuous",
      region_predicates = list(function(A, H) rep(TRUE, length(A))),
      policy_pieces = list(function(A, H) A),
      inverse_map_pieces = list(function(A_star, H) A_star),
      inverse_deriv_pieces = list(function(A_star, H) rep(1, length(A_star))),
      name = "identity"
    )
  } else {
    mtp_id <- mtp_discrete(
      map_fun = function(A, H) A,
      support = discrete_support,
      name = "identity"
    )
  }

  repeat_policy_over_time(mtp_id, tau = tau, name = "identity_sequence")
}

.make_subgroup_matrix <- function(ds, subgroup_funs) {
  if (is.function(subgroup_funs)) {
    subgroup_funs <- list(group = subgroup_funs)
  }
  if (!is.list(subgroup_funs) || !all(vapply(subgroup_funs, is.function, logical(1)))) {
    stop("`subgroup_funs` must be a function or named list of functions.")
  }
  out <- lapply(subgroup_funs, function(f) as.numeric(f(ds$df)))
  out <- as.data.frame(out, check.names = FALSE)
  if (is.null(names(out))) {
    names(out) <- paste0("group", seq_len(ncol(out)))
  }
  out
}

.make_scalar_H_fun <- function(k_provider, t) {
  function(A_vec, H_df) {
    # current scalar path only uses observed-data K_obs(t)
    # returned as a constant-in-evaluation-point function for the observed update
    k_provider$K_obs(t)
  }
}

.make_subgroup_H_obs <- function(k_provider, subgroup_mat, t) {
  pA <- colMeans(subgroup_mat)
  pA <- pmax(pA, 1e-8)
  sweep(subgroup_mat, 2, pA, "/") * k_provider$K_obs(t)
}

# helpers for fluctuation in subgroup LMTP --------------------------------

# Evaluate subgroup indicators on an arbitrary data.frame H_df
.make_subgroup_matrix_from_df <- function(H_df, subgroup_funs) {
  if (is.function(subgroup_funs)) {
    subgroup_funs <- list(group = subgroup_funs)
  }
  if (!is.list(subgroup_funs) || !all(vapply(subgroup_funs, is.function, logical(1)))) {
    stop("`subgroup_funs` must be a function or named list of functions.")
  }

  out <- lapply(subgroup_funs, function(f) {
    val <- as.numeric(f(H_df))
    if (length(val) != nrow(H_df)) {
      stop("Each subgroup function must return a vector of length nrow(H_df).")
    }
    if (!all(val %in% c(0, 1))) {
      stop("Each subgroup function must return 0/1 indicators.")
    }
    val
  })

  out <- as.data.frame(out, check.names = FALSE)
  if (is.null(names(out))) {
    names(out) <- paste0("group", seq_len(ncol(out)))
  }
  out
}

.make_subgroup_H_fun <- function(subgroup_funs, pA, K_prev_obs, r_t_fun) {
  force(subgroup_funs)
  force(pA)
  force(K_prev_obs)
  force(r_t_fun)

  function(A_vec, H_df) {
    if (length(A_vec) != nrow(H_df)) {
      stop("`A_vec` must have length nrow(H_df).")
    }

    G_eval <- .make_subgroup_matrix_from_df(H_df, subgroup_funs)
    G_eval_scaled <- sweep(G_eval, 2, pA, "/")

    ratio_eval <- K_prev_obs * r_t_fun(A_vec, H_df)

    # multiply each subgroup column by the scalar ratio vector
    sweep(G_eval_scaled, 1, ratio_eval, "*")
  }
}


#' Clip Probabilities Away from 0-1 by a Small Amount
.clip_probability <- function(x, clip_probability = 1e-6) {
  pmin(pmax(x, clip_probability), 1 - clip_probability)
}

#' Lift Probabilities away from 0 -- Especially for Density Estimates
.truncate_positive <- function(x, truncate_density = 1e-12) {
  pmax(x, truncate_density)
}


#' Truncate to an Interval -- Especially for Density Ratios
.truncate_interval <- function(x, truncate = c(1e-6, Inf)) {
  pmin(pmax(x, truncate[1]), truncate[2])
}


