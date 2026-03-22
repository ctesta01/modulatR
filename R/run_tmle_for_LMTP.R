# core runner: scalar E[Y^d] -------------------------------------------------

#' Run a scalar TMLE for an LMTP mean
#'
#' @param ds An `LMTPData` object.
#' @param policy_seq An `LMTPPolicySequence`.
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#' @param learners_Q Learner spec passed to `q_trainer()`.
#' @param fml_Q Optional outcome-regression formula/formulas.
#' @param learners_Q_extra_args Optional extra args for Q learners.
#' @param alpha Confidence level parameter.
#'
#' @return An `LMTPFit`.
#' @export
run_tmle_for_LMTP <- function(ds,
                              policy_seq,
                              nuisance_factory,
                              fluctuation,
                              learners_Q,
                              fml_Q = NULL,
                              learners_Q_extra_args = NULL,
                              alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) stop("`ds` must inherit from `LMTPData`.")
  if (!inherits(policy_seq, "LMTPPolicySequence")) stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) stop("`fluctuation` must inherit from `LMTPFluctuationSubmodel`.")

  nuisance_factory$train(ds)
  kprov <- nuisance_factory$k_provider(ds)
  tau <- ds$tau()

  q_fit_factory <- nuisance_factory$q_trainer(
    learners_Q = learners_Q,
    fml_Q = fml_Q,
    learners_Q_extra_args = learners_Q_extra_args
  )

  Q_init <- vector("list", tau)
  Q_star <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)

  target_next <- ds$Y()

  for (t in rev(seq_len(tau))) {
    Q0_t <- q_fit_factory(ds, t = t, pseudo_outcome_vec = target_next)
    Qvec_t <- Q0_t(ds$A(t), ds$H(t))

    H_obs_t <- kprov$K_obs(t)
    H_fun_t <- .make_scalar_H_fun(kprov, t)

    up <- fluctuation$fit_update(
      Q0_fun_t = Q0_t,
      Qvec = Qvec_t,
      target_vec = target_next,
      H_obs = H_obs_t,
      H_fun_t = H_fun_t,
      t = t
    )

    Q_init[[t]] <- Q0_t
    Q_star[[t]] <- up$wrap
    eps[[t]] <- up$epsilon
    intercept[[t]] <- up$intercept

    if (t > 1L) {
      A_t_star <- policy_seq$apply_t(t, ds$A(t), ds$H(t))
      target_next <- Q_star[[t]](A_t_star, ds$H(t))
    }
  }

  A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))
  psi <- mean(Q_star[[1]](A1_star, ds$H(1)))

  ic <- rep(0, ds$n)
  for (t in seq_len(tau)) {
    A_t_star <- policy_seq$apply_t(t, ds$A(t), ds$H(t))
    Q_t_star <- Q_star[[t]](A_t_star, ds$H(t))

    if (t < tau) {
      A_next_star <- policy_seq$apply_t(t + 1L, ds$A(t + 1L), ds$H(t + 1L))
      Q_next_star <- Q_star[[t + 1L]](A_next_star, ds$H(t + 1L))
    } else {
      Q_next_star <- ds$Y()
    }

    ic <- ic + kprov$K_obs(t) * (Q_next_star - Q_t_star)
  }
  ic <- ic + Q_star[[1]](A1_star, ds$H(1)) - psi

  var_hat <- stats::var(ic) / ds$n
  se_hat <- sqrt(var_hat)
  ci_hat <- .wald_ci(psi, se_hat, alpha = alpha)

  LMTPFit$new(
    estimate = psi,
    var = var_hat,
    se = se_hat,
    ci = ci_hat,
    alpha = alpha,
    parameter = "E[Y^d]",
    estimator = "LMTP TMLE",
    n = ds$n,
    ic = ic,
    eif = ic,
    Q_init = Q_init,
    Q_star = Q_star,
    omega = kprov$K_obs_all,
    eps = eps,
    intercept = intercept
  )
}

# subgroup runner -------------------------------------------------------------

#' Run a simultaneous subgroup TMLE for LMTP subgroup means
#'
#' @param ds An `LMTPData` object.
#' @param policy_seq An `LMTPPolicySequence`.
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#' @param subgroup_funs A function or named list of functions taking `ds$df`
#'   and returning subgroup indicators.
#' @param learners_Q Learner spec passed to `q_trainer()`.
#' @param fml_Q Optional formula/formulas for `Q_t`.
#' @param learners_Q_extra_args Optional extra args for Q learners.
#' @param alpha Confidence level parameter.
#'
#' @return An `LMTPFit` with vector estimate and EIF matrix.
#' @export
run_subgroup_LMTP_TMLE <- function(ds,
                                   policy_seq,
                                   nuisance_factory,
                                   fluctuation,
                                   subgroup_funs,
                                   learners_Q,
                                   fml_Q = NULL,
                                   learners_Q_extra_args = NULL,
                                   alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) stop("`ds` must inherit from `LMTPData`.")
  if (!inherits(policy_seq, "LMTPPolicySequence")) stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) stop("`fluctuation` must inherit from `LMTPFluctuationSubmodel`.")

  subgroup_mat <- .make_subgroup_matrix(ds, subgroup_funs)
  subgroup_names <- colnames(subgroup_mat)
  n_groups <- ncol(subgroup_mat)

  nuisance_factory$train(ds)
  kprov <- nuisance_factory$k_provider(ds)
  tau <- ds$tau()

  q_fit_factory <- nuisance_factory$q_trainer(
    learners_Q = learners_Q,
    fml_Q = fml_Q,
    learners_Q_extra_args = learners_Q_extra_args
  )

  Q_init <- vector("list", tau)
  Q_star <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)

  target_next <- ds$Y()

  for (t in rev(seq_len(tau))) {
    Q0_t <- q_fit_factory(ds, t = t, pseudo_outcome_vec = target_next)
    Qvec_t <- Q0_t(ds$A(t), ds$H(t))

    H_obs_t <- .make_subgroup_H_obs(kprov, subgroup_mat, t)

    # current implementation keeps observed-point subgroup H for fitting.
    # evaluation-time subgroup H is left open for later refinement when
    # subgroup LMTP clever covariates are fully formalized.
    H_fun_t <- function(A_vec, H_df) H_obs_t

    up <- fluctuation$fit_update(
      Q0_fun_t = Q0_t,
      Qvec = Qvec_t,
      target_vec = target_next,
      H_obs = H_obs_t,
      H_fun_t = H_fun_t,
      t = t
    )

    Q_init[[t]] <- Q0_t
    Q_star[[t]] <- up$wrap
    eps[[t]] <- up$epsilon
    intercept[[t]] <- up$intercept

    if (t > 1L) {
      A_t_star <- policy_seq$apply_t(t, ds$A(t), ds$H(t))
      target_next <- Q_star[[t]](A_t_star, ds$H(t))
    }
  }

  A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))
  Q1_star <- Q_star[[1]](A1_star, ds$H(1))

  pA <- pmax(colMeans(subgroup_mat), 1e-8)
  psi <- colMeans(sweep(subgroup_mat, 2, pA, "/") * Q1_star)

  eif <- matrix(0, nrow = ds$n, ncol = n_groups)
  colnames(eif) <- subgroup_names

  base_ic <- rep(0, ds$n)
  for (t in seq_len(tau)) {
    A_t_star <- policy_seq$apply_t(t, ds$A(t), ds$H(t))
    Q_t_star <- Q_star[[t]](A_t_star, ds$H(t))

    if (t < tau) {
      A_next_star <- policy_seq$apply_t(t + 1L, ds$A(t + 1L), ds$H(t + 1L))
      Q_next_star <- Q_star[[t + 1L]](A_next_star, ds$H(t + 1L))
    } else {
      Q_next_star <- ds$Y()
    }

    base_ic <- base_ic + kprov$K_obs(t) * (Q_next_star - Q_t_star)
  }

  for (j in seq_len(n_groups)) {
    Sj <- subgroup_mat[[j]] / pA[j]
    eif[, j] <- Sj * base_ic + Sj * Q1_star - psi[j]
  }

  var_hat <- apply(eif, 2, stats::var) / ds$n
  se_hat <- sqrt(var_hat)
  ci_hat <- t(vapply(
    seq_len(n_groups),
    function(j) .wald_ci(psi[j], se_hat[j], alpha = alpha),
    numeric(2)
  ))
  rownames(ci_hat) <- subgroup_names

  LMTPFit$new(
    estimate = psi,
    var = var_hat,
    se = se_hat,
    ci = ci_hat,
    alpha = alpha,
    parameter = "Subgroup E[Y^d]",
    estimator = "Subgroup LMTP TMLE",
    n = ds$n,
    ic = eif,
    eif = eif,
    Q_init = Q_init,
    Q_star = Q_star,
    omega = kprov$K_obs_all,
    eps = eps,
    intercept = intercept,
    subgroup_names = subgroup_names
  )
}

# TL_Task wrappers ------------------------------------------------------------

#' Construct a TL_Task for scalar LMTP TMLE
#'
#' @param ds An `LMTPData` object.
#' @param policy_seq An `LMTPPolicySequence`.
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#' @param learners_Q Learner spec for the Q recursion.
#' @param fml_Q Optional formulas for Q.
#' @param learners_Q_extra_args Optional extra args for Q learners.
#'
#' @return A `TL_Task`.
#' @export
make_lmtp_tl_task <- function(ds,
                              policy_seq,
                              nuisance_factory,
                              fluctuation,
                              learners_Q,
                              fml_Q = NULL,
                              learners_Q_extra_args = NULL) {

  task <- TL_Task$new(
    parameter = list(name = "E[Y^d]"),
    data = ds,
    nuisance_factory = nuisance_factory,
    fluctuation = fluctuation,
    targeting_step = function() NULL,
    run = function(alpha = 0.05) {
      self$results <- run_tmle_for_LMTP(
        ds = self$data,
        policy_seq = self$policy_seq,
        nuisance_factory = self$nuisances,
        fluctuation = self$fluctuation,
        learners_Q = self$learners_Q,
        fml_Q = self$fml_Q,
        learners_Q_extra_args = self$learners_Q_extra_args,
        alpha = alpha
      )
      self$results
    }
  )

  task$policy_seq <- policy_seq
  task$learners_Q <- learners_Q
  task$fml_Q <- fml_Q
  task$learners_Q_extra_args <- learners_Q_extra_args

  task
}

#' Construct a TL_Task for subgroup LMTP TMLE
#'
#' @param ds An `LMTPData` object.
#' @param policy_seq An `LMTPPolicySequence`.
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#' @param subgroup_funs Function or named list of subgroup indicator functions.
#' @param learners_Q Learner spec for the Q recursion.
#' @param fml_Q Optional formulas for Q.
#' @param learners_Q_extra_args Optional extra args for Q learners.
#'
#' @return A `TL_Task`.
#' @export
make_subgroup_lmtp_tl_task <- function(ds,
                                       policy_seq,
                                       nuisance_factory,
                                       fluctuation,
                                       subgroup_funs,
                                       learners_Q,
                                       fml_Q = NULL,
                                       learners_Q_extra_args = NULL) {

  task <- TL_Task$new(
    parameter = list(name = "Subgroup E[Y^d]"),
    data = ds,
    nuisance_factory = nuisance_factory,
    fluctuation = fluctuation,
    targeting_step = function() NULL,
    run = function(alpha = 0.05) {
      self$results <- run_subgroup_LMTP_TMLE(
        ds = self$data,
        policy_seq = self$policy_seq,
        nuisance_factory = self$nuisances,
        fluctuation = self$fluctuation,
        subgroup_funs = self$subgroup_funs,
        learners_Q = self$learners_Q,
        fml_Q = self$fml_Q,
        learners_Q_extra_args = self$learners_Q_extra_args,
        alpha = alpha
      )
      self$results
    }
  )

  task$policy_seq <- policy_seq
  task$subgroup_funs <- subgroup_funs
  task$learners_Q <- learners_Q
  task$fml_Q <- fml_Q
  task$learners_Q_extra_args <- learners_Q_extra_args

  task
}

# delta-method contrasts ------------------------------------------------------

#' Delta-method contrast for two LMTP fits
#'
#' @param fit1 An `LMTPFit`.
#' @param fit0 An `LMTPFit`.
#' @param transform One of `"difference"`, `"rr"`, `"or"`.
#' @param alpha Confidence level parameter.
#'
#' @return An `LMTPFit`.
#' @export
lmtp_delta_contrast <- function(fit1,
                                fit0,
                                transform = c("difference", "rr", "or"),
                                alpha = 0.05) {
  transform <- match.arg(transform)

  if (!inherits(fit1, "LMTPFit") || !inherits(fit0, "LMTPFit")) {
    stop("`fit1` and `fit0` must inherit from `LMTPFit`.")
  }
  if (length(fit1$estimate) != 1L || length(fit0$estimate) != 1L) {
    stop("`lmtp_delta_contrast()` currently expects scalar fits.")
  }
  if (length(fit1$eif) != length(fit0$eif)) {
    stop("Fits must have EIFs of the same length.")
  }

  psi1 <- fit1$estimate
  psi0 <- fit0$estimate
  ic1 <- fit1$eif
  ic0 <- fit0$eif

  if (transform == "difference") {
    est <- psi1 - psi0
    ic <- ic1 - ic0
    param_name <- "E[Y^d] - E[Y]"
  } else if (transform == "rr") {
    est <- psi1 / psi0
    ic <- (1 / psi0) * ic1 - (psi1 / psi0^2) * ic0
    param_name <- "E[Y^d] / E[Y]"
  } else {
    odds <- function(p) p / (1 - p)
    est <- odds(psi1) / odds(psi0)
    g1 <- 1 / (psi1 * (1 - psi1))
    g0 <- -1 / (psi0 * (1 - psi0))
    log_or_ic <- g1 * ic1 + g0 * ic0
    ic <- est * log_or_ic
    param_name <- "OR(E[Y^d], E[Y])"
  }

  var_hat <- stats::var(ic) / length(ic)
  se_hat <- sqrt(var_hat)
  ci_hat <- .wald_ci(est, se_hat, alpha = alpha)

  LMTPFit$new(
    estimate = est,
    var = var_hat,
    se = se_hat,
    ci = ci_hat,
    alpha = alpha,
    parameter = param_name,
    estimator = paste("Delta-method", transform),
    n = length(ic),
    ic = ic,
    eif = ic
  )
}
