# core runner: E[Y^d] -------------------------------------------------

#' Estimate a scalar TMLE for an LMTP mean
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
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }
  if (!inherits(policy_seq, "LMTPPolicySequence")) {
    stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
  }
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) {
    stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
  }
  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) {
    stop("`fluctuation` must inherit from `LMTPFluctuationSubmodel`.")
  }

  policy_seq$validate_against_data(ds)
  nuisance_factory$train(ds)
  kprov <- nuisance_factory$k_provider(ds)

  tau <- ds$tau()

  q_fit_factory <- nuisance_factory$q_trainer(
    learners_Q = learners_Q,
    fml_Q = fml_Q,
    learners_Q_extra_args = learners_Q_extra_args
  )

  # start of backward recursions --------------------------------------------

  Q_init <- vector("list", tau)
  Q_star <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)

  # -----------------------------
  # Track 1: untargeted recursion
  # -----------------------------
  target_next <- ds$Y()

  for (t in rev(seq_len(tau))) {
    A_t <- ds$A(t)
    H_t <- ds$H(t)
    A_t_star <- policy_seq$apply_t(t, A_t, H_t)

    Q0_t <- q_fit_factory(ds, t = t, pseudo_outcome_vec = target_next)

    Q_init[[t]] <- Q0_t

    if (t > 1L) {
      target_next <- Q_init[[t]](A_t_star, H_t)
    }
  }

  # ----------------------------------------
  # Track 2: targeted recursion for TMLE psi
  # ----------------------------------------
  target_next_star <- ds$Y()

  for (t in rev(seq_len(tau))) {
    A_t <- ds$A(t)
    H_t <- ds$H(t)
    A_t_star <- policy_seq$apply_t(t, A_t, H_t)

    Q0_t <- Q_init[[t]]
    Qvec_t <- Q0_t(A_t, H_t)

    H_obs_t <- kprov$K_obs(t)

    K_prev_obs <- if (t == 1L) rep(1, ds$n) else kprov$K_obs(t - 1L)

    r_t_fun <- nuisance_factory$r_list[[t]]
    if (is.null(r_t_fun) || !is.function(r_t_fun)) {
      stop("`r_t_fun` is NULL or not a function at time t = ", t, ".")
    }

    H_fun_t <- local({
      K_prev_obs_local <- K_prev_obs
      r_t_fun_local <- r_t_fun
      function(A_vec, H_df) {
        K_prev_obs_local * r_t_fun_local(A_vec, H_df)
      }
    })

    up <- fluctuation$fit_update(
      Q0_fun_t = Q0_t,
      Qvec = Qvec_t,
      target_vec = target_next_star,
      H_obs = H_obs_t,
      H_fun_t = H_fun_t,
      t = t
    )

    Q_star[[t]] <- up$wrap
    eps[[t]] <- up$epsilon
    intercept[[t]] <- up$intercept

    if (t > 1L) {
      target_next_star <- Q_star[[t]](A_t_star, H_t)
    }
  }

  # Plug-in estimate
  A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))

  psi_tmle <- mean(Q_star[[1]](A1_star, ds$H(1)))
  psi_plugin_init <- mean(Q_init[[1]](A1_star, ds$H(1)))

  # EIF matching the parameterization:
  #
  # for t < tau:
  #   K_t { Q_{t+1}(A_{t+1}^d, H_{t+1}) - Q_t(A_t, H_t) }
  #
  # for t = tau:
  #   K_tau { Y - Q_tau(A_tau, H_tau) }
  #
  # plus plug-in term:
  #   Q_1(A_1^d, H_1) - psi
  ic <- rep(0, ds$n)

  for (t in seq_len(tau)) {
    K_t <- kprov$K_obs(t)

    if (t < tau) {
      A_next_star <- policy_seq$apply_t(t + 1L, ds$A(t + 1L), ds$H(t + 1L))

      Q_next_init_shift <- Q_init[[t + 1L]](A_next_star, ds$H(t + 1L))
      Q_t_init_obs <- Q_init[[t]](ds$A(t), ds$H(t))

      ic <- ic + K_t * (Q_next_init_shift - Q_t_init_obs)
    } else {
      Q_t_init_obs <- Q_init[[t]](ds$A(t), ds$H(t))
      ic <- ic + K_t * (ds$Y() - Q_t_init_obs)
    }
  }

  ic <- ic + Q_init[[1]](A1_star, ds$H(1)) - psi_plugin_init

  var_hat <- stats::var(ic) / ds$n
  se_hat <- sqrt(var_hat)
  ci_hat <- .wald_ci(psi_tmle, se_hat, alpha = alpha)

  LMTPFit$new(
    estimate = psi_tmle,
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

# run_subgroup_tmle_for_LMTP -------------------------------------------------------------

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
run_subgroup_tmle_for_LMTP <- function(ds,
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

  policy_seq$validate_against_data(ds)

  if (isTRUE(nuisance_factory$cross_fit)) {
    stop(
      "`run_subgroup_LMTP_TMLE()` currently requires `cross_fit = FALSE` ",
      "until the subgroup fluctuation path is fully implemented for cross-fitting."
    )
  }

  subgroup_mat <- .make_subgroup_matrix(ds, subgroup_funs)

  if (nrow(subgroup_mat) != ds$n) {
    stop("`subgroup_funs` must return vectors of length `ds$n`.")
  }
  if (!all(vapply(subgroup_mat, function(x) all(is.finite(x)), logical(1)))) {
    stop("Subgroup indicators must be finite.")
  }
  if (!all(vapply(subgroup_mat, function(x) all(x %in% c(0, 1)), logical(1)))) {
    stop("Each subgroup function must return a 0/1 indicator.")
  }

  pA <- colMeans(subgroup_mat)
  if (any(pA <= 0)) {
    bad <- colnames(subgroup_mat)[pA <= 0]
    stop("The following subgroup(s) have zero empirical prevalence: ",
         paste(bad, collapse = ", "))
  }

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

  # sequential backwards regressions -----------------------------------------

  Q_init <- vector("list", tau)
  Q_star <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)

  # --------------------------------
  # Track 1: untargeted Q recursion
  # --------------------------------
  target_next <- ds$Y()

  for (t in rev(seq_len(tau))) {
    H_t <- ds$H(t)
    A_t <- ds$A(t)
    A_t_star <- policy_seq$apply_t(t, A_t, H_t)

    Q0_t <- q_fit_factory(ds, t = t, pseudo_outcome_vec = target_next)
    Q_init[[t]] <- Q0_t

    if (t > 1L) {
      target_next <- Q_init[[t]](A_t_star, H_t)
    }
  }

  # --------------------------------
  # Track 2: targeted Q recursion
  # --------------------------------
  target_next_star <- ds$Y()

  for (t in rev(seq_len(tau))) {
    H_t <- ds$H(t)
    A_t <- ds$A(t)
    A_t_star <- policy_seq$apply_t(t, A_t, H_t)

    Q0_t <- Q_init[[t]]
    Qvec_t <- Q0_t(A_t, H_t)

    K_prev_obs <- if (t == 1L) rep(1, ds$n) else kprov$K_obs(t - 1L)

    r_t_fun <- nuisance_factory$r_list[[t]]
    if (is.null(r_t_fun) || !is.function(r_t_fun)) {
      stop("`r_t_fun` is NULL or not a function at time t = ", t, ".")
    }

    H_fun_t <- .make_subgroup_H_fun(
      subgroup_funs = subgroup_funs,
      pA = pA,
      K_prev_obs = K_prev_obs,
      r_t_fun = r_t_fun
    )

    H_obs_t <- H_fun_t(A_t, H_t)

    up <- fluctuation$fit_update(
      Q0_fun_t = Q0_t,
      Qvec = Qvec_t,
      target_vec = target_next_star,
      H_obs = H_obs_t,
      H_fun_t = H_fun_t,
      t = t
    )

    Q_star[[t]] <- up$wrap
    eps[[t]] <- up$epsilon
    intercept[[t]] <- up$intercept

    if (t > 1L) {
      target_next_star <- Q_star[[t]](A_t_star, H_t)
    }
  }

  A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))

  Q1_star_shift <- Q_star[[1]](A1_star, ds$H(1))
  Q1_init_shift <- Q_init[[1]](A1_star, ds$H(1))

  weights_subgroup <- sweep(subgroup_mat, 2, pA, "/")

  psi_tmle <- colMeans(weights_subgroup * Q1_star_shift)
  psi_plugin_init <- colMeans(weights_subgroup * Q1_init_shift)

  # eif construction --------------------------------------------------------

  pA <- pmax(colMeans(subgroup_mat), 1e-8)

  eif <- matrix(0, nrow = ds$n, ncol = n_groups)
  colnames(eif) <- subgroup_names

  ic <- rep(0, ds$n)
  ic_star <- rep(0, ds$n)

  for (t in seq_len(tau)) {
    K_t <- kprov$K_obs(t)

    if (t < tau) {
      H_next <- ds$H(t + 1L)
      A_next_star <- policy_seq$apply_t(t + 1L, ds$A(t + 1L), H_next)

      Q_next_init_shift <- Q_init[[t + 1L]](A_next_star, H_next)
      Q_t_init_obs <- Q_init[[t]](ds$A(t), ds$H(t))
      ic <- ic + K_t * (Q_next_init_shift - Q_t_init_obs)

      Q_next_star_shift <- Q_star[[t + 1L]](A_next_star, H_next)
      Q_t_star_obs <- Q_star[[t]](ds$A(t), ds$H(t))
      ic_star <- ic_star + K_t * (Q_next_star_shift - Q_t_star_obs)
    } else {
      Q_t_init_obs <- Q_init[[t]](ds$A(t), ds$H(t))
      ic <- ic + K_t * (ds$Y() - Q_t_init_obs)

      Q_t_star_obs <- Q_star[[t]](ds$A(t), ds$H(t))
      ic_star <- ic_star + K_t * (ds$Y() - Q_t_star_obs)
    }
  }

  eif <- matrix(0, nrow = ds$n, ncol = n_groups)
  eif_star <- matrix(0, nrow = ds$n, ncol = n_groups)
  colnames(eif) <- subgroup_names
  colnames(eif_star) <- subgroup_names

  for (g in seq_len(n_groups)) {
    wg <- subgroup_mat[, g] / pA[g]

    eif[, g] <- wg * ic + wg * Q1_init_shift - psi_plugin_init[g]
    eif_star[, g] <- wg * ic_star + wg * Q1_star_shift - psi_tmle[g]
  }

  var_hat <- apply(eif_star, 2, stats::var) / ds$n
  se_hat <- sqrt(var_hat)

  ci_hat <- rbind(
    lower = psi_tmle - stats::qnorm(1 - alpha / 2) * se_hat,
    upper = psi_tmle + stats::qnorm(1 - alpha / 2) * se_hat
  )
  colnames(ci_hat) <- subgroup_names


  LMTPFit$new(
    estimate = psi_tmle,
    var = var_hat,
    se = se_hat,
    ci = ci_hat,
    alpha = alpha,
    parameter = "Subgroup E[Y^d]",
    estimator = "Subgroup LMTP TMLE",
    n = ds$n,
    ic = eif_star,
    eif = eif_star,
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
