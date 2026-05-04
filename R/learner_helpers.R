#' Repeat a learner constructor over all time points
#'
#' @param learner_fun A learner function taking one `LMTPData` argument.
#' @param tau Number of treatment time points.
#'
#' @return A list of learner functions of length `tau`.
#' @export
repeat_lmtp_learner <- function(learner_fun, tau) {
  if (!is.function(learner_fun)) {
    stop("`learner_fun` must be a function.")
  }
  if (!is.numeric(tau) || length(tau) != 1L || tau < 1L) {
    stop("`tau` must be a positive scalar.")
  }

  rep(list(learner_fun), as.integer(tau))
}

#' Construct GLM learners for all LMTP outcome regressions
#'
#' @param tau Number of treatment time points.
#' @param family GLM family.
#' @param formula Optional formula or list of formulas. If `NULL`, a formula
#'   using all columns in `(A_t, H_t)` is constructed.
#' @param exclude_id Logical; whether to exclude the id column from regressors.
#' @param ... Additional arguments passed to `stats::glm()`.
#'
#' @return A list of learner functions.
#' @export
make_glm_m_learners <- function(tau,
                                family = stats::gaussian(),
                                formula = NULL,
                                exclude_id = TRUE,
                                ...) {
  repeat_lmtp_learner(
    make_glm_m_learner(
      family = family,
      formula = formula,
      exclude_id = exclude_id,
      ...
    ),
    tau = tau
  )
}

#' Construct a GLM learner for one LMTP sequential regression m_t
#'
#' @param family GLM family.
#' @param formula Optional formula or list of formulas.
#' @param exclude_id Logical; whether to exclude the id column from regressors.
#' @param ... Additional arguments passed to `stats::glm()`.
#'
#' @export
make_glm_m_learner <- function(family = stats::gaussian(),
                               formula = NULL,
                               exclude_id = TRUE,
                               ...) {
  glm_args <- list(...)

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    t <- data$metadata$t
    pseudo_outcome_col <- data$metadata$pseudo_outcome_col
    policy_seq <- data$metadata$policy_seq

    if (is.null(t)) stop("`data$metadata$t` is missing.")
    if (is.null(pseudo_outcome_col)) {
      stop("`data$metadata$pseudo_outcome_col` is missing.")
    }
    if (is.null(policy_seq)) {
      stop("`data$metadata$policy_seq` is missing.")
    }

    A_name <- data$A_cols[[t]]
    H_t <- data$H(t)
    A_t <- data$A(t)
    A_t_d <- policy_seq$apply_t(t, A_t, H_t)

    dat <- data$AH(t)
    dat[[pseudo_outcome_col]] <- data$df[[pseudo_outcome_col]]

    fml <- private_formula_for_time(
      formula = formula,
      t = t,
      outcome_col = pseudo_outcome_col,
      regressors = colnames(data$AH(t)),
      id_col = data$id_col,
      exclude_id = exclude_id
    )

    fit <- do.call(
      stats::glm,
      c(
        list(formula = fml, family = family, data = dat),
        glm_args
      )
    )

    nd_obs <- H_t
    nd_obs[[A_name]] <- A_t

    nd_d <- H_t
    nd_d[[A_name]] <- A_t_d

    list(
      fit = fit,
      formula = fml,
      m_obs = as.numeric(stats::predict(fit, newdata = nd_obs, type = "response")),
      m_d = as.numeric(stats::predict(fit, newdata = nd_d, type = "response"))
    )
  }
}

#' Construct nadir density learners for all LMTP treatment mechanisms
#'
#' @param tau Number of treatment time points.
#' @param learner A nadir learner function.
#' @param formula Optional formula or list of formulas.
#' @param exclude_id Logical; whether to exclude the id column from regressors.
#' @param ... Additional arguments passed to the nadir learner.
#'
#' @return A list of learner functions.
#' @export
make_nadir_density_g_learners <- function(tau,
                                          learner = nadir::lnr_glm_density,
                                          formula = NULL,
                                          exclude_id = TRUE,
                                          ...) {
  repeat_lmtp_learner(
    make_nadir_density_g_learner(
      learner = learner,
      formula = formula,
      exclude_id = exclude_id,
      ...
    ),
    tau = tau
  )
}

#' Construct a nadir density learner for one LMTP treatment mechanism g_t
#'
#' @param learner A nadir learner function.
#' @param formula Optional formula or list of formulas.
#' @param exclude_id Logical; whether to exclude the id column from regressors.
#' @param ... Additional arguments passed to the nadir learner.
#'
#' @return A learner function taking one `LMTPData` object.
#' @export
make_nadir_density_g_learner <- function(learner = nadir::lnr_glm_density,
                                         formula = NULL,
                                         exclude_id = TRUE,
                                         ...) {
  learner_args <- list(...)

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    t <- data$metadata$t
    if (is.null(t)) {
      stop("`data$metadata$t` is missing.")
    }

    A_name <- data$A_cols[[t]]
    H_t <- data$H(t)
    dat <- data$AH(t)

    fml <- private_formula_for_time(
      formula = formula,
      t = t,
      outcome_col = A_name,
      regressors = colnames(H_t),
      id_col = data$id_col,
      exclude_id = exclude_id
    )

    fit <- do.call(
      learner,
      c(
        list(
          formula = fml,
          data = dat
        ),
        learner_args
      )
    )

    predict_density_fun <- function(A_vec, H_df) {
      if (length(A_vec) != nrow(H_df)) {
        stop("`A_vec` must have length `nrow(H_df)`.")
      }

      newdata <- H_df
      newdata[[A_name]] <- A_vec

      if (is.function(fit)) {
        return(as.numeric(fit(newdata)))
      }

      as.numeric(stats::predict(fit, newdata = newdata))
    }

    list(
      fit = fit,
      formula = fml,
      predict_density = predict_density_fun
    )
  }
}

make_density_ratio_learner <- function(g_learner,
                                       truncate_density = 1e-12) {
  force(g_learner)

  function(data) {
    t <- data$metadata$t
    policy_seq <- data$metadata$policy_seq

    if (is.null(t)) stop("`data$metadata$t` is missing.")
    if (is.null(policy_seq)) {
      stop("`data$metadata$policy_seq` is missing.")
    }

    A_t <- data$A(t)
    H_t <- data$H(t)
    A_t_d <- policy_seq$apply_t(t, A_t, H_t)

    g_fit <- g_learner(data)

    if (is.null(g_fit$predict_density) ||
        !is.function(g_fit$predict_density)) {
      stop("`g_learner` must return a list with `predict_density`.")
    }

    g_fun <- g_fit$predict_density

    g_obs <- .truncate_positive(
      g_fun(A_t, H_t),
      truncate = truncate_density
    )

    g_d_obs <- .truncate_positive(
      policy_seq$gd_t(
        t = t,
        A_vec = A_t,
        H_df = H_t,
        density_fun = g_fun
      ),
      truncate = truncate_density
    )

    g_at_d <- .truncate_positive(
      g_fun(A_t_d, H_t),
      truncate = truncate_density
    )

    g_d_at_d <- .truncate_positive(
      policy_seq$gd_t(
        t = t,
        A_vec = A_t_d,
        H_df = H_t,
        density_fun = g_fun
      ),
      truncate = truncate_density
    )

    list(
      fit = g_fit,
      r_obs = as.numeric(g_d_obs / g_obs),
      r_d = as.numeric(g_d_at_d / g_at_d)
    )
  }
}

private_formula_for_time <- function(formula,
                                     t,
                                     outcome_col,
                                     regressors,
                                     id_col = NULL,
                                     exclude_id = TRUE) {
  if (!is.null(formula)) {
    if (inherits(formula, "formula")) {
      return(formula)
    }

    if (is.list(formula) && length(formula) >= t && inherits(formula[[t]], "formula")) {
      return(formula[[t]])
    }

    stop("`formula` must be NULL, a formula, or a list of formulas.")
  }

  rhs <- regressors

  if (isTRUE(exclude_id) && !is.null(id_col)) {
    rhs <- setdiff(rhs, id_col)
  }

  if (length(rhs) == 0L) {
    return(stats::as.formula(paste(outcome_col, "~ 1")))
  }

  stats::as.formula(
    paste(outcome_col, "~", paste(rhs, collapse = " + "))
  )
}
