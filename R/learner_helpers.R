# LMTP learner helpers -----------------------------------------------------

repeat_lmtp_learner <- function(learner_fun, tau) {
  if (!is.function(learner_fun)) {
    stop("`learner_fun` must be a function.")
  }
  if (!is.numeric(tau) || length(tau) != 1L || tau < 1L) {
    stop("`tau` must be a positive scalar.")
  }

  rep(list(learner_fun), as.integer(tau))
}

#' Construct GLM m learners for all time points
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

#' Construct one GLM m learner
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

    if (is.null(t)) {
      stop("`data$metadata$t` is missing.")
    }
    if (is.null(pseudo_outcome_col)) {
      stop("`data$metadata$pseudo_outcome_col` is missing.")
    }

    fit_df <- data$AH(t)
    fit_df[[pseudo_outcome_col]] <- data$df[[pseudo_outcome_col]]

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
        list(
          formula = fml,
          family = family,
          data = fit_df
        ),
        glm_args
      )
    )

    function(newdata) {
      if (!inherits(newdata, "LMTPData")) {
        stop("`newdata` must inherit from `LMTPData`.")
      }

      as.numeric(stats::predict(
        fit,
        newdata = newdata$AH(t),
        type = "response"
      ))
    }
  }
}

# Nadir m learner ----------------------------------------------------------

#' Construct one nadir m learner
#'
#' @export
make_nadir_m_learner <- function(learner = nadir::lnr_ranger,
                                 formula = NULL,
                                 exclude_id = TRUE,
                                 ...) {
  learner_args <- list(...)

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    t <- data$metadata$t
    pseudo_outcome_col <- data$metadata$pseudo_outcome_col

    if (is.null(t)) {
      stop("`data$metadata$t` is missing.")
    }
    if (is.null(pseudo_outcome_col)) {
      stop("`data$metadata$pseudo_outcome_col` is missing.")
    }

    fit_df <- data$AH(t)
    fit_df[[pseudo_outcome_col]] <- data$df[[pseudo_outcome_col]]

    fml <- private_formula_for_time(
      formula = formula,
      t = t,
      outcome_col = pseudo_outcome_col,
      regressors = colnames(data$AH(t)),
      id_col = data$id_col,
      exclude_id = exclude_id
    )

    fit <- do.call(
      learner,
      c(
        list(
          formula = fml,
          data = fit_df
        ),
        learner_args
      )
    )

    function(newdata) {
      if (!inherits(newdata, "LMTPData")) {
        stop("`newdata` must inherit from `LMTPData`.")
      }

      pred_df <- newdata$AH(t)

      if (is.function(fit)) {
        return(as.numeric(fit(pred_df)))
      }

      as.numeric(stats::predict(fit, newdata = pred_df))
    }
  }
}


#' Construct one nadir density g learner
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
    if (is.null(t)) stop("`data$metadata$t` is missing.")

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
        list(formula = fml, data = dat),
        learner_args
      )
    )

    function(newdata) {
      if (!inherits(newdata, "LMTPData")) {
        stop("`newdata` must inherit from `LMTPData`.")
      }

      pred_df <- newdata$AH(t)

      if (is.function(fit)) {
        return(as.numeric(fit(pred_df)))
      }

      as.numeric(stats::predict(fit, newdata = pred_df))
    }
  }
}

#' Construct nadir density g learners for all time points
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


# Binary lambda learners ---------------------------------------------------

#' Construct one GLM lambda learner for ratio classification
#'
#' @export
make_glm_lambda_learner <- function(family = stats::binomial(),
                                    formula = NULL,
                                    exclude_id = TRUE,
                                    ...) {
  glm_args <- list(...)

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    t <- data$metadata$t
    lambda_col <- data$metadata$lambda_col

    if (is.null(t)) {
      stop("`data$metadata$t` is missing.")
    }
    if (is.null(lambda_col)) {
      stop("`data$metadata$lambda_col` is missing.")
    }

    fit_df <- data$AH(t)
    fit_df[[lambda_col]] <- data$df[[lambda_col]]

    fml <- private_formula_for_time(
      formula = formula,
      t = t,
      outcome_col = lambda_col,
      regressors = colnames(data$AH(t)),
      id_col = data$id_col,
      exclude_id = exclude_id
    )

    fit <- do.call(
      stats::glm,
      c(
        list(
          formula = fml,
          family = family,
          data = fit_df
        ),
        glm_args
      )
    )

    function(newdata) {
      if (!inherits(newdata, "LMTPData")) {
        stop("`newdata` must inherit from `LMTPData`.")
      }

      as.numeric(stats::predict(
        fit,
        newdata = newdata$AH(t),
        type = "response"
      ))
    }
  }
}


#' Construct GLM lambda learners for all time points
#'
#' @export
make_glm_lambda_learners <- function(tau,
                                     family = stats::binomial(),
                                     formula = NULL,
                                     exclude_id = TRUE,
                                     ...) {
  repeat_lmtp_learner(
    make_glm_lambda_learner(
      family = family,
      formula = formula,
      exclude_id = exclude_id,
      ...
    ),
    tau = tau
  )
}
