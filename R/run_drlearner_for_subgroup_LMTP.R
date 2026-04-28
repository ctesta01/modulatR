# ============================================================
# DR-learner summaries for subgroup-specific LMTP means
# ============================================================

# Helper: make scalar DR pseudo-outcome
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


# Helper: fit subgroup means of a pseudo-outcome
#
# For mutually exclusive/exhaustive subgroup indicators and no intercept,
# OLS is equivalent to estimating subgroup means.
fit_subgroup_dr_means <- function(ds,
                                  pseudo_outcome,
                                  subgroup_funs) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }
  if (!is.numeric(pseudo_outcome) || length(pseudo_outcome) != ds$n) {
    stop("`pseudo_outcome` must be numeric of length ds$n.")
  }

  subgroup_mat <- .make_subgroup_matrix(ds, subgroup_funs)
  subgroup_names <- colnames(subgroup_mat)

  if (!all(vapply(subgroup_mat, function(x) all(x %in% c(0, 1)), logical(1)))) {
    stop("Subgroup functions must return 0/1 indicators.")
  }

  # make syntactically valid names for regression
  safe_names <- make.names(subgroup_names, unique = TRUE)
  subgroup_mat_safe <- subgroup_mat
  colnames(subgroup_mat_safe) <- safe_names

  dat <- data.frame(pseudo_outcome = pseudo_outcome, subgroup_mat_safe, check.names = FALSE)

  # no intercept: one coefficient per subgroup
  fml <- stats::as.formula(
    paste("pseudo_outcome ~ -1 +", paste(safe_names, collapse = " + "))
  )

  fit <- stats::lm(fml, data = dat)

  pA <- colMeans(subgroup_mat)
  n_g <- colSums(subgroup_mat)

  est <- numeric(length(subgroup_names))
  se <- numeric(length(subgroup_names))

  for (g in seq_along(subgroup_names)) {
    idx <- subgroup_mat[, g] == 1
    est[g] <- mean(pseudo_outcome[idx])
    se[g] <- stats::sd(pseudo_outcome[idx]) / sqrt(sum(idx))
  }

  names(est) <- subgroup_names
  names(se) <- subgroup_names

  list(
    subgroup_names = subgroup_names,
    safe_names = safe_names,
    subgroup_mat = subgroup_mat,
    prevalence = pA,
    n_g = n_g,
    lm_fit = fit,
    estimate = est,
    se = se,
    ci_lower = est - stats::qnorm(0.975) * se,
    ci_upper = est + stats::qnorm(0.975) * se
  )
}


# Main function:
#   E[Y^d | subgroup], E[Y | subgroup], E[Y^d - Y | subgroup]
#
# using scalar TMLE fits + DR pseudo-outcomes + subgroup-level regression
run_drlearner_tmle_lmtp_subgroup_summary <- function(ds,
                                                     policy_seq,
                                                     nuisance_factory,
                                                     fluctuation,
                                                     subgroup_funs,
                                                     learners_Q,
                                                     control_policy_seq = NULL,
                                                     fml_Q = NULL,
                                                     learners_Q_extra_args = NULL,
                                                     alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) stop("`ds` must inherit from `LMTPData`.")
  if (!inherits(policy_seq, "LMTPPolicySequence")) stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) stop("`fluctuation` must inherit from `LMTPFluctuationSubmodel`.")

  if (is.null(control_policy_seq)) {
    control_policy_seq <- identity_policy(
      A_type = "continuous",
      tau = ds$tau()
    )
  }

  # --------------------------------------------
  # Scalar TMLE fit for E[Y^d]
  # --------------------------------------------
  fit_d <- run_tmle_for_LMTP(
    ds = ds,
    policy_seq = policy_seq,
    nuisance_factory = nuisance_factory,
    fluctuation = fluctuation,
    learners_Q = learners_Q,
    fml_Q = fml_Q,
    learners_Q_extra_args = learners_Q_extra_args,
    alpha = alpha
  )

  # --------------------------------------------
  # Scalar TMLE fit for E[Y]
  # --------------------------------------------
  fit_0 <- run_tmle_for_LMTP(
    ds = ds,
    policy_seq = control_policy_seq,
    nuisance_factory = nuisance_factory,
    fluctuation = fluctuation,
    learners_Q = learners_Q,
    fml_Q = fml_Q,
    learners_Q_extra_args = learners_Q_extra_args,
    alpha = alpha
  )

  # --------------------------------------------
  # DR pseudo-outcomes
  # --------------------------------------------
  phi_d <- make_lmtp_dr_pseudooutcome(fit_d)
  phi_0 <- make_lmtp_dr_pseudooutcome(fit_0)
  phi_diff <- phi_d - phi_0

  # --------------------------------------------
  # subgroup-specific DR means
  # --------------------------------------------
  out_d <- fit_subgroup_dr_means(
    ds = ds,
    pseudo_outcome = phi_d,
    subgroup_funs = subgroup_funs
  )

  out_0 <- fit_subgroup_dr_means(
    ds = ds,
    pseudo_outcome = phi_0,
    subgroup_funs = subgroup_funs
  )

  out_diff <- fit_subgroup_dr_means(
    ds = ds,
    pseudo_outcome = phi_diff,
    subgroup_funs = subgroup_funs
  )

  subgroup_names <- out_d$subgroup_names

  summary_table <- data.frame(
    subgroup = subgroup_names,

    EYd_est = as.numeric(out_d$estimate[subgroup_names]),
    EYd_se  = as.numeric(out_d$se[subgroup_names]),

    EY_est  = as.numeric(out_0$estimate[subgroup_names]),
    EY_se   = as.numeric(out_0$se[subgroup_names]),

    contrast_est = as.numeric(out_diff$estimate[subgroup_names]),
    contrast_se  = as.numeric(out_diff$se[subgroup_names]),

    prevalence = as.numeric(out_d$prevalence[subgroup_names]),
    n_g = as.numeric(out_d$n_g[subgroup_names]),

    row.names = NULL
  )

  obj <- list(
    fit_d = fit_d,
    fit_0 = fit_0,
    phi_d = phi_d,
    phi_0 = phi_0,
    phi_diff = phi_diff,
    dr_EYd = out_d,
    dr_EY = out_0,
    dr_contrast = out_diff,
    summary_table = summary_table
  )
  class(obj) <- "LMTPDRSubgroupSummary"
  obj
}


print.LMTPDRSubgroupSummary <- function(x, ...) {
  cat("LMTPDRSubgroupSummary\n\n")

  tab <- x$summary_table

  for (i in seq_len(nrow(tab))) {
    sg <- tab$subgroup[i]

    cat(sg, "\n", sep = "")
    cat("  E[Y^d | ", sg, "]      = ",
        formatC(tab$EYd_est[i], digits = 4, format = "f"),
        "  (SE = ",
        formatC(tab$EYd_se[i], digits = 4, format = "f"),
        ")\n", sep = "")

    cat("  E[Y | ", sg, "]        = ",
        formatC(tab$EY_est[i], digits = 4, format = "f"),
        "  (SE = ",
        formatC(tab$EY_se[i], digits = 4, format = "f"),
        ")\n", sep = "")

    cat("  E[Y^d - Y | ", sg, "]  = ",
        formatC(tab$contrast_est[i], digits = 4, format = "f"),
        "  (SE = ",
        formatC(tab$contrast_se[i], digits = 4, format = "f"),
        ")\n\n", sep = "")
  }

  invisible(x)
}


fit_lmtp_dr_meta_learner <- function(ds,
                                     pseudo_outcome,
                                     subgroup_funs = NULL,
                                     meta_formula = NULL) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }
  if (!is.numeric(pseudo_outcome) || length(pseudo_outcome) != ds$n) {
    stop("`pseudo_outcome` must be numeric of length ds$n.")
  }

  if (!is.null(subgroup_funs)) {
    X <- .make_subgroup_matrix(ds, subgroup_funs)
    original_names <- colnames(X)
    safe_names <- make.names(original_names, unique = TRUE)
    colnames(X) <- safe_names
  } else {
    X <- ds$W()
    if (is.null(X)) {
      X <- data.frame()
    }
    original_names <- colnames(X)
    safe_names <- colnames(X)
  }

  meta_dat <- data.frame(pseudo_outcome = pseudo_outcome, X, check.names = FALSE)

  if (is.null(meta_formula)) {
    if (!is.null(subgroup_funs)) {
      meta_formula <- stats::as.formula(
        paste("pseudo_outcome ~ -1 +", paste(safe_names, collapse = " + "))
      )
    } else if (ncol(X) == 0L) {
      meta_formula <- pseudo_outcome ~ 1
    } else {
      meta_formula <- stats::as.formula("pseudo_outcome ~ .")
    }
  }

  meta_fit <- stats::lm(meta_formula, data = meta_dat)
  fitted_vals <- as.numeric(stats::predict(meta_fit, newdata = meta_dat))

  list(
    meta_fit = meta_fit,
    meta_formula = meta_formula,
    meta_data = meta_dat,
    fitted_values = fitted_vals,
    design_matrix = X,
    original_names = original_names,
    safe_names = safe_names
  )
}
