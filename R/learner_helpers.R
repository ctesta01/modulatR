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


# Binary treatment-mass GLM learner ----

#' Construct a GLM g learner for a binary treatment
#'
#' @description
#' Constructs a treatment-mechanism learner for binary discrete treatments.
#' The learner fits a binomial generalized linear model for
#' `P(A_t = 1 | H_t)` (more generally, for one of the two treatment levels)
#' and returns a prediction function that evaluates the conditional treatment
#' mass
#'
#' `g_t(A_t | H_t)`.
#'
#' Thus, unlike an ordinary binomial prediction function, the returned
#' prediction function does not always return the probability of the
#' "success" category. It returns the probability corresponding to the
#' treatment value contained in each row of `newdata`.
#'
#' This makes the learner compatible with `LMTPNuisanceFactory`, which must
#' evaluate the treatment mass at both observed and intervention-modified
#' treatment values.
#'
#' @param family A binomial GLM family. Defaults to logistic regression via
#'   `stats::binomial()`.
#' @param formula Optional treatment-model formula. May be a single formula
#'   or a list of time-specific formulas. If `NULL`, the treatment is regressed
#'   on all variables in `H_t`, excluding the ID column by default.
#' @param exclude_id Logical; if `TRUE`, exclude the ID column from an
#'   automatically constructed formula.
#' @param ... Additional arguments passed to `stats::glm()`.
#'
#' @return A learner function with contract
#'   `learner(data) -> function(newdata)`, where the prediction function
#'   returns `g_t(A_t | H_t)`.
#'
#' @export
make_glm_g_learner <- function(family = stats::binomial(),
                               formula = NULL,
                               exclude_id = TRUE,
                               ...) {
  glm_args <- list(...)

  if (!inherits(family, "family") ||
      !family$family %in% c("binomial", "quasibinomial")) {
    stop(
      "`make_glm_g_learner()` requires a binomial or quasibinomial family."
    )
  }

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    t <- data$metadata$t

    if (is.null(t)) {
      stop("`data$metadata$t` is missing.")
    }

    A_name <- data$A_cols[[t]]
    A_train <- data$A(t)

    if (anyNA(A_train)) {
      stop(
        "`make_glm_g_learner()` does not support missing treatment values."
      )
    }

    # Determine the two observed treatment levels.
    if (is.factor(A_train)) {
      support <- levels(droplevels(A_train))
    } else {
      support <- sort(unique(A_train))
    }

    if (length(support) != 2L) {
      stop(
        "`make_glm_g_learner()` requires a binary treatment. ",
        "At time ", t, ", found ", length(support),
        " observed treatment level(s)."
      )
    }

    support_key <- as.character(support)

    # The choice of which level is coded as "success" is arbitrary as long
    # as fitting and prediction use the same convention.
    success_key <- support_key[[2L]]

    binary_outcome_col <- "..g_binary_outcome"

    fit_df <- data$H(t)
    fit_df[[binary_outcome_col]] <-
      as.numeric(as.character(A_train) == success_key)

    # Construct the treatment-model formula using the existing package helper.
    #
    # We first build/read the formula with A_t as its nominal outcome so that
    # user-supplied formulas such as
    #
    #   A1 ~ W1 + W2
    #
    # have the expected public interface. We then replace the left-hand side
    # with the internal 0/1 representation used by glm().
    fml <- private_formula_for_time(
      formula = formula,
      t = t,
      outcome_col = A_name,
      regressors = colnames(data$H(t)),
      id_col = data$id_col,
      exclude_id = exclude_id
    )

    if (length(fml) != 3L) {
      stop(
        "`formula` must be a two-sided formula such as `A1 ~ W1 + W2`."
      )
    }

    rhs_vars <- all.vars(stats::delete.response(stats::terms(fml)))

    if (A_name %in% rhs_vars) {
      stop(
        "The current treatment variable `", A_name,
        "` cannot appear on the right-hand side of its own g model."
      )
    }

    fml[[2L]] <- as.name(binary_outcome_col)

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

      A_new <- newdata$A(t)
      A_new_key <- as.character(A_new)

      if (anyNA(A_new)) {
        stop(
          "Cannot evaluate the treatment mass at missing treatment values."
        )
      }

      outside_support <- !A_new_key %in% support_key

      if (any(outside_support)) {
        bad_values <- unique(A_new_key[outside_support])

        stop(
          "Encountered treatment value(s) outside the fitted binary support: ",
          paste(bad_values, collapse = ", "),
          ". Fitted support was: ",
          paste(support_key, collapse = ", "),
          "."
        )
      }

      p_success <- as.numeric(
        stats::predict(
          fit,
          newdata = newdata$H(t),
          type = "response"
        )
      )

      is_success <- A_new_key == success_key

      ifelse(
        is_success,
        p_success,
        1 - p_success
      )
    }
  }
}
