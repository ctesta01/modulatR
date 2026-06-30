# DR-learner for conditional LMTP means ------------------------------------

#' Construct a DR pseudo-outcome from a scalar LMTP fit
#'
#' @param fit A scalar `LMTPFit`.
#'
#' @return Numeric vector of length `fit$n`.
#' @export
make_lmtp_dr_pseudooutcome <- function(fit) {
  if (!inherits(fit, "LMTPFit")) {
    stop("`fit` must inherit from `LMTPFit`.")
  }

  if (length(fit$estimate) != 1L) {
    stop("`fit` must be scalar.")
  }

  if (is.null(fit$eif)) {
    stop("`fit$eif` is NULL.")
  }

  as.numeric(fit$estimate + fit$eif)
}


# Meta-learner helpers ------------------------------------------------------

#' Construct a simple linear-model meta-learner
#'
#' @description
#' A meta-learner takes a data.frame, an outcome column, and covariate columns,
#' then returns a prediction function. This default implementation uses `lm()`.
#'
#' @param formula Optional formula. If `NULL`, uses `outcome_col ~ .`.
#' @param ... Additional arguments passed to `stats::lm()`.
#'
#' @return A meta-learner function.
#' @export
make_lm_metalearner <- function(formula = NULL, ...) {
  lm_args <- list(...)

  function(data, outcome_col, covariate_cols) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.")
    }
    if (!outcome_col %in% colnames(data)) {
      stop("`outcome_col` is not present in `data`.")
    }
    if (!all(covariate_cols %in% colnames(data))) {
      missing <- setdiff(covariate_cols, colnames(data))
      stop("Missing covariate columns: ", paste(missing, collapse = ", "))
    }

    fit_data <- data[, c(outcome_col, covariate_cols), drop = FALSE]

    if (is.null(formula)) {
      if (length(covariate_cols) == 0L) {
        fml <- stats::as.formula(paste(outcome_col, "~ 1"))
      } else {
        fml <- stats::as.formula(paste(outcome_col, "~ ."))
      }
    } else {
      fml <- formula
    }

    fit <- do.call(
      stats::lm,
      c(
        list(formula = fml, data = fit_data),
        lm_args
      )
    )

    list(
      fit = fit,
      predict = function(newdata) {
        if (!is.data.frame(newdata)) {
          stop("`newdata` must be a data.frame.")
        }

        as.numeric(stats::predict(fit, newdata = newdata))
      }
    )
  }
}


#' Construct a nadir meta-learner
#'
#' @description
#' This is useful for continuous conditioning covariates, where the second-stage
#' regression can be nonlinear or nonparametric.
#'
#' @param learner A nadir learner, for example `nadir::lnr_ranger`.
#' @param formula Optional formula.
#' @param ... Additional arguments passed to `learner`.
#'
#' @return A meta-learner function.
#' @export
make_nadir_metalearner <- function(learner,
                                   formula = NULL,
                                   ...) {
  learner_args <- list(...)

  function(data, outcome_col, covariate_cols) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.")
    }
    if (!outcome_col %in% colnames(data)) {
      stop("`outcome_col` is not present in `data`.")
    }
    if (!all(covariate_cols %in% colnames(data))) {
      missing <- setdiff(covariate_cols, colnames(data))
      stop("Missing covariate columns: ", paste(missing, collapse = ", "))
    }

    fit_data <- data[, c(outcome_col, covariate_cols), drop = FALSE]

    if (is.null(formula)) {
      if (length(covariate_cols) == 0L) {
        fml <- stats::as.formula(paste(outcome_col, "~ 1"))
      } else {
        fml <- stats::as.formula(paste(outcome_col, "~ ."))
      }
    } else {
      fml <- formula
    }

    fit <- do.call(
      learner,
      c(
        list(formula = fml, data = fit_data),
        learner_args
      )
    )

    list(
      fit = fit,
      predict = function(newdata) {
        if (!is.data.frame(newdata)) {
          stop("`newdata` must be a data.frame.")
        }

        if (is.function(fit)) {
          return(as.numeric(fit(newdata)))
        }

        as.numeric(stats::predict(fit, newdata = newdata))
      }
    )
  }
}


# Conditioning data helpers ------------------------------------------------

.make_conditioning_data <- function(ds,
                                    conditioning_cols = NULL,
                                    conditioning_fun = NULL) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!is.null(conditioning_fun)) {
    V <- conditioning_fun(ds$df)
  } else if (!is.null(conditioning_cols)) {
    V <- ds$df[, conditioning_cols, drop = FALSE]
  } else {
    V <- ds$W()
    if (is.null(V)) {
      V <- data.frame(intercept_only = rep(1, ds$n))
    }
  }

  V <- as.data.frame(V, check.names = FALSE)

  if (nrow(V) != ds$n) {
    stop("Conditioning data must have `ds$n` rows.")
  }

  if (is.null(colnames(V)) || any(colnames(V) == "")) {
    colnames(V) <- paste0("V", seq_len(ncol(V)))
  }

  V
}


#' Fit a meta-learner to an LMTP DR pseudo-outcome
#'
#' @param ds An `LMTPData` object.
#' @param pseudo_outcome Numeric vector of length `ds$n`.
#' @param metalearner A function with signature
#'   `metalearner(data, outcome_col, covariate_cols)`.
#' @param conditioning_cols Optional column names in `ds$df`.
#' @param conditioning_fun Optional function taking `ds$df` and returning a
#'   data.frame of conditioning covariates.
#'
#' @return A list containing the fitted meta-learner and fitted values.
#' @export
fit_lmtp_dr_metalearner <- function(ds,
                                    pseudo_outcome,
                                    metalearner = make_lm_metalearner(),
                                    conditioning_cols = NULL,
                                    conditioning_fun = NULL) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!is.numeric(pseudo_outcome) || length(pseudo_outcome) != ds$n) {
    stop("`pseudo_outcome` must be numeric of length `ds$n`.")
  }

  if (!is.function(metalearner)) {
    stop("`metalearner` must be a function.")
  }

  V <- .make_conditioning_data(
    ds = ds,
    conditioning_cols = conditioning_cols,
    conditioning_fun = conditioning_fun
  )

  outcome_col <- "..dr_pseudooutcome"
  meta_data <- data.frame(
    ..dr_pseudooutcome = pseudo_outcome,
    V,
    check.names = FALSE
  )

  covariate_cols <- colnames(V)

  meta_fit <- metalearner(
    data = meta_data,
    outcome_col = outcome_col,
    covariate_cols = covariate_cols
  )

  if (!is.list(meta_fit) || !is.function(meta_fit$predict)) {
    stop("`metalearner` must return a list with a `$predict` function.")
  }

  fitted_values <- as.numeric(meta_fit$predict(meta_data))

  if (length(fitted_values) != ds$n) {
    stop("The meta-learner prediction function returned the wrong length.")
  }

  list(
    meta_fit = meta_fit,
    meta_data = meta_data,
    conditioning_data = V,
    covariate_cols = covariate_cols,
    fitted_values = fitted_values
  )
}


# Nuisance factory cloning -------------------------------------------------
#
# The identity/no-shift fit needs its own nuisance factory because the policy
# lives in the factory. This helper copies the learner specifications and
# numerical truncation settings, but swaps the policy.

clone_lmtp_nuisance_factory_with_policy <- function(nuisance_factory,
                                                    policy_seq) {
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) {
    stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  }
  if (!inherits(policy_seq, "LMTPPolicySequence")) {
    stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
  }

  ratio_args <- if (!is.null(nuisance_factory$g_learners)) {
    list(g_learners = nuisance_factory$g_learners)
  } else if (!is.null(nuisance_factory$lambda_learners)) {
    list(lambda_learners = nuisance_factory$lambda_learners)
  } else {
    stop("`nuisance_factory` must contain either `g_learners` or `lambda_learners`.")
  }

  do.call(
    LMTPNuisanceFactory$new,
    c(
      list(
        policy_seq = policy_seq,
        m_learners = nuisance_factory$m_learners,
        truncate_density = nuisance_factory$truncate_density,
        truncate_ratio = nuisance_factory$truncate_ratio,
        clip_lambda_probability = nuisance_factory$clip_lambda_probability
      ),
      ratio_args
    )
  )
}


# Subgroup summaries from pseudo-outcomes ----------------------------------

fit_subgroup_dr_means <- function(ds,
                                  pseudo_outcome,
                                  subgroup_funs,
                                  alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!is.numeric(pseudo_outcome) || length(pseudo_outcome) != ds$n) {
    stop("`pseudo_outcome` must be numeric of length `ds$n`.")
  }

  subgroup_mat <- .make_subgroup_matrix(ds, subgroup_funs)
  subgroup_mat <- as.data.frame(subgroup_mat, check.names = FALSE)

  if (nrow(subgroup_mat) != ds$n) {
    stop("`subgroup_funs` must return vectors of length `ds$n`.")
  }

  subgroup_mat[] <- lapply(subgroup_mat, as.numeric)

  if (!all(vapply(subgroup_mat, function(x) all(x %in% c(0, 1)), logical(1)))) {
    stop("Each subgroup function must return a 0/1 indicator.")
  }

  subgroup_names <- colnames(subgroup_mat)
  if (is.null(subgroup_names) || any(subgroup_names == "")) {
    subgroup_names <- paste0("subgroup", seq_len(ncol(subgroup_mat)))
    colnames(subgroup_mat) <- subgroup_names
  }

  prevalence <- colMeans(subgroup_mat)
  n_g <- colSums(subgroup_mat)

  if (any(n_g == 0)) {
    bad <- subgroup_names[n_g == 0]
    stop("The following subgroup(s) have no observations: ", paste(bad, collapse = ", "))
  }

  est <- numeric(length(subgroup_names))
  se <- numeric(length(subgroup_names))

  for (g in seq_along(subgroup_names)) {
    idx <- subgroup_mat[[g]] == 1
    est[g] <- mean(pseudo_outcome[idx])
    se[g] <- stats::sd(pseudo_outcome[idx]) / sqrt(sum(idx))
  }

  names(est) <- subgroup_names
  names(se) <- subgroup_names

  z <- stats::qnorm(1 - alpha / 2)

  list(
    subgroup_names = subgroup_names,
    subgroup_mat = subgroup_mat,
    prevalence = prevalence,
    n_g = n_g,
    estimate = est,
    se = se,
    ci_lower = est - z * se,
    ci_upper = est + z * se
  )
}


# Main DR-learner runner ---------------------------------------------------

#' Run a DR-learner for conditional LMTP means and contrasts
#'
#' @description
#' Fits scalar LMTP TMLEs, constructs doubly robust pseudo-outcomes, and then
#' regresses those pseudo-outcomes on user-specified conditioning variables.
#'
#' Unlike `run_subgroup_tmle_for_LMTP()`, this function is not restricted to a
#' finite vector of subgroup indicators. Its main use case is estimating
#' conditional LMTP functions over continuous or multivariate conditioning
#' covariates.
#'
#' @param ds An `LMTPData` object.
#' @param nuisance_factory An `LMTPNuisanceFactory` for the shifted policy.
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#' @param metalearner A second-stage regression learner with signature
#'   `metalearner(data, outcome_col, covariate_cols)`.
#' @param conditioning_cols Optional column names in `ds$df` defining `V`.
#' @param conditioning_fun Optional function taking `ds$df` and returning a
#'   data.frame of conditioning covariates.
#' @param control_nuisance_factory Optional nuisance factory for the control or
#'   identity policy.
#' @param control_policy_seq Optional policy used to construct the control
#'   nuisance factory when `control_nuisance_factory = NULL`.
#' @param A_type Treatment type passed to `identity_policy()` if no control
#'   policy is supplied.
#' @param subgroup_funs Optional subgroup functions. If supplied, subgroup means
#'   of the DR pseudo-outcomes are also returned.
#' @param alpha Wald interval level for subgroup summaries and scalar fits.
#'
#' @return An object of class `LMTPDRLearner`.
#' @export
run_drlearner_for_subgroup_LMTP <- function(ds,
                                            nuisance_factory,
                                            fluctuation = LMTPFluctuationSubmodel$new(),
                                            metalearner = make_lm_metalearner(),
                                            conditioning_cols = NULL,
                                            conditioning_fun = NULL,
                                            control_nuisance_factory = NULL,
                                            control_policy_seq = NULL,
                                            A_type = "continuous",
                                            subgroup_funs = NULL,
                                            alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) {
    stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  }

  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) {
    stop("`fluctuation` must inherit from `LMTPFluctuationSubmodel`.")
  }

  if (!is.function(metalearner)) {
    stop("`metalearner` must be a function.")
  }

  if (is.null(control_nuisance_factory)) {
    if (is.null(control_policy_seq)) {
      control_policy_seq <- identity_policy(
        A_type = A_type,
        tau = ds$tau()
      )
    }

    control_nuisance_factory <- clone_lmtp_nuisance_factory_with_policy(
      nuisance_factory = nuisance_factory,
      policy_seq = control_policy_seq
    )
  }

  if (!inherits(control_nuisance_factory, "LMTPNuisanceFactory")) {
    stop("`control_nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  }

  # Scalar LMTP fits -------------------------------------------------------

  fit_d <- run_tmle_for_LMTP(
    ds = ds,
    nuisance_factory = nuisance_factory,
    fluctuation = fluctuation,
    alpha = alpha
  )

  fit_0 <- run_tmle_for_LMTP(
    ds = ds,
    nuisance_factory = control_nuisance_factory,
    fluctuation = fluctuation,
    alpha = alpha
  )

  # DR pseudo-outcomes -----------------------------------------------------

  phi_d <- make_lmtp_dr_pseudooutcome(fit_d)
  phi_0 <- make_lmtp_dr_pseudooutcome(fit_0)
  phi_diff <- phi_d - phi_0

  # Second-stage meta-learning --------------------------------------------

  meta_d <- fit_lmtp_dr_metalearner(
    ds = ds,
    pseudo_outcome = phi_d,
    metalearner = metalearner,
    conditioning_cols = conditioning_cols,
    conditioning_fun = conditioning_fun
  )

  meta_0 <- fit_lmtp_dr_metalearner(
    ds = ds,
    pseudo_outcome = phi_0,
    metalearner = metalearner,
    conditioning_cols = conditioning_cols,
    conditioning_fun = conditioning_fun
  )

  meta_diff <- fit_lmtp_dr_metalearner(
    ds = ds,
    pseudo_outcome = phi_diff,
    metalearner = metalearner,
    conditioning_cols = conditioning_cols,
    conditioning_fun = conditioning_fun
  )

  # Optional subgroup summaries -------------------------------------------

  subgroup_summary <- NULL

  if (!is.null(subgroup_funs)) {
    out_d <- fit_subgroup_dr_means(
      ds = ds,
      pseudo_outcome = phi_d,
      subgroup_funs = subgroup_funs,
      alpha = alpha
    )

    out_0 <- fit_subgroup_dr_means(
      ds = ds,
      pseudo_outcome = phi_0,
      subgroup_funs = subgroup_funs,
      alpha = alpha
    )

    out_diff <- fit_subgroup_dr_means(
      ds = ds,
      pseudo_outcome = phi_diff,
      subgroup_funs = subgroup_funs,
      alpha = alpha
    )

    subgroup_names <- out_d$subgroup_names

    subgroup_summary <- data.frame(
      subgroup = subgroup_names,

      EYd_est = as.numeric(out_d$estimate[subgroup_names]),
      EYd_se = as.numeric(out_d$se[subgroup_names]),

      EY_est = as.numeric(out_0$estimate[subgroup_names]),
      EY_se = as.numeric(out_0$se[subgroup_names]),

      contrast_est = as.numeric(out_diff$estimate[subgroup_names]),
      contrast_se = as.numeric(out_diff$se[subgroup_names]),

      prevalence = as.numeric(out_d$prevalence[subgroup_names]),
      n_g = as.numeric(out_d$n_g[subgroup_names]),

      row.names = NULL
    )
  }

  out <- list(
    fit_d = fit_d,
    fit_0 = fit_0,

    phi_d = phi_d,
    phi_0 = phi_0,
    phi_diff = phi_diff,

    meta_EYd = meta_d,
    meta_EY = meta_0,
    meta_contrast = meta_diff,

    subgroup_summary = subgroup_summary,

    conditioning_cols = conditioning_cols,
    conditioning_fun = conditioning_fun,
    alpha = alpha
  )

  class(out) <- "LMTPDRLearner"
  out
}


# Prediction method --------------------------------------------------------

#' Predict from an LMTP DR-learner
#'
#' @param object An `LMTPDRLearner`.
#' @param newdata A data.frame containing the conditioning covariates.
#' @param type One of `"EYd"`, `"EY"`, or `"contrast"`.
#' @param ... Ignored.
#'
#' @return Numeric predictions.
#' @export
predict.LMTPDRLearner <- function(object,
                                  newdata,
                                  type = c("EYd", "EY", "contrast"),
                                  ...) {
  type <- match.arg(type)

  if (!inherits(object, "LMTPDRLearner")) {
    stop("`object` must inherit from `LMTPDRLearner`.")
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  fit <- switch(
    type,
    EYd = object$meta_EYd$meta_fit,
    EY = object$meta_EY$meta_fit,
    contrast = object$meta_contrast$meta_fit
  )

  as.numeric(fit$predict(newdata))
}


# Print method -------------------------------------------------------------

#' @export
print.LMTPDRLearner <- function(x, ...) {
  cat("LMTPDRLearner\n")
  cat("  scalar shifted estimate E[Y^d]: ", formatC(x$fit_d$estimate, digits = 4, format = "f"), "\n", sep = "")
  cat("  scalar control estimate E[Y]:   ", formatC(x$fit_0$estimate, digits = 4, format = "f"), "\n", sep = "")
  cat("  scalar contrast:                ", formatC(x$fit_d$estimate - x$fit_0$estimate, digits = 4, format = "f"), "\n", sep = "")

  if (!is.null(x$subgroup_summary)) {
    cat("\nSubgroup DR summaries:\n")
    print(x$subgroup_summary, row.names = FALSE)
  }

  invisible(x)
}
