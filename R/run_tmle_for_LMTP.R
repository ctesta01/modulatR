
# this file handles the core runners for estimating E[Y^d]  and E[Y^d | V = v]

#' Estimate a scalar TMLE for an LMTP mean
#'
#' @description
#' Vector-based TMLE runner for
#' \eqn{\psi = E[Y(\bar A^d)]}.
#'
#' All nuisance specifications live in `nuisance_factory`.
#'
#' The plug-in estimate uses targeted regressions `m_star`, while the estimated
#' EIF used for variance uses untargeted regressions `m_init`.
#'
#' @param ds An `LMTPData` object.
#' @param nuisance_factory An `LMTPNuisanceFactory` object.
#' @param fluctuation An `LMTPFluctuationSubmodel` object.
#' @param alpha Wald interval level.
#'
#' @return An `LMTPFit`.
#' @export
run_tmle_for_LMTP <- function(ds,
                              nuisance_factory,
                              fluctuation = LMTPFluctuationSubmodel$new(),
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

  policy_seq <- nuisance_factory$policy_seq
  if (!inherits(policy_seq, "LMTPPolicySequence")) {
    stop("`nuisance_factory$policy_seq` must inherit from `LMTPPolicySequence`.")
  }

  policy_seq$validate_against_data(ds)

  tau <- ds$tau()
  n <- ds$n

  #
  # 1. Fit ratio nuisance vectors.
  #
  # Stores:
  #   r_obs[[t]]     = r_t(A_t, H_t)
  #   r_d[[t]]       = r_t(A_t^d, H_t)
  #   omega_obs[[t]] = prod_{s <= t} r_s(A_s, H_s)
  #

  nuisance_factory$train_ratios(ds)

  #
  # 2. Fit initial, untargeted sequential regressions.
  #
  # Stores:
  #   m_init_obs[[t]] = m_t(A_t, H_t)
  #   m_init_d[[t]]   = m_t(A_t^d, H_t)
  #
  # The pseudo-outcome recursion is:
  #   pseudo_tau = Y
  #   pseudo_{t-1} = m_t(A_t^d, H_t)
  #

  m_init_obs <- vector("list", tau)
  m_init_d <- vector("list", tau)

  pseudo_outcome_init <- ds$Y()

  for (t in rev(seq_len(tau))) {
    m_t <- nuisance_factory$train_m_t(
      ds = ds,
      t = t,
      pseudo_outcome = pseudo_outcome_init
    )

    m_init_obs[[t]] <- m_t$m_obs
    m_init_d[[t]] <- m_t$m_d

    if (t > 1L) {
      pseudo_outcome_init <- m_t$m_d
    }
  }

  #
  # 3. Targeted backward recursion.
  #
  # The fluctuation object handles the submodel. The runner supplies vectors:
  #
  #   observed clever covariate:
  #     H_obs_t = omega_t
  #
  #   policy-evaluation clever covariate:
  #     H_d_t = omega_{t-1} * r_t(A_t^d, H_t)
  #
  # The plug-in recursion uses targeted policy-value predictions m_t^d,*.
  #

  m_star_obs <- vector("list", tau)
  m_star_d <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)
  fluctuation_fits <- vector("list", tau)

  pseudo_outcome_star <- ds$Y()

  for (t in rev(seq_len(tau))) {
    H_obs_t <- nuisance_factory$omega(t)
    H_d_t <- nuisance_factory$omega_prev(t, n = n) * nuisance_factory$r_d_preds[[t]]

    update_t <- fluctuation$fit_update(
      m_obs = m_init_obs[[t]],
      target = pseudo_outcome_star,
      H_obs = H_obs_t,
      m_d = m_init_d[[t]],
      H_d = H_d_t,
      t = t
    )

    m_star_obs[[t]] <- update_t$m_obs_star
    m_star_d[[t]] <- update_t$m_d_star
    eps[[t]] <- update_t$epsilon
    intercept[[t]] <- update_t$intercept
    fluctuation_fits[[t]] <- update_t$fit

    if (t > 1L) {
      pseudo_outcome_star <- update_t$m_d_star
    }
  }

  #
  # 4. TMLE plug-in estimate.
  #
  # Uses targeted m_1 evaluated at the modified treatment value.
  #

  psi_hat <- mean(m_star_d[[1L]])

  #
  # 5. Estimated EIF / IC.
  #
  # Crucial: use untargeted nuisance regressions m_init, not m_star.
  #
  #   D_hat(O_i)
  #     = sum_t omega_t(O_i)
  #         { next_untargeted_t(O_i) - m_t(A_t,H_t)(O_i) }
  #       + m_1(A_1^d,H_1)(O_i) - psi_hat
  #
  # where:
  #   next_untargeted_t = Y for t = tau
  #   next_untargeted_t = m_{t+1}(A_{t+1}^d,H_{t+1}) otherwise.
  #

  eif <- rep(0, n)

  for (t in seq_len(tau)) {
    next_term_untargeted <- if (t == tau) {
      ds$Y()
    } else {
      m_init_d[[t + 1L]]
    }

    eif <- eif +
      nuisance_factory$omega(t) *
      (next_term_untargeted - m_init_obs[[t]])
  }

  eif <- eif + m_init_d[[1L]] - psi_hat

  var_hat <- stats::var(eif) / n
  se_hat <- sqrt(var_hat)
  ci <- .wald_ci(psi_hat, se_hat, alpha = alpha)

  LMTPFit$new(
    estimate = psi_hat,
    var = var_hat,
    se = se_hat,
    ci = ci,
    alpha = alpha,
    parameter = "E[Y^d]",
    estimator = "TMLE",
    n = n,
    ic = eif,
    eif = eif,

    m_init = list(
      m_obs = m_init_obs,
      m_d = m_init_d
    ),
    m_star = list(
      m_obs = m_star_obs,
      m_d = m_star_d
    ),
    omega = nuisance_factory$omega_preds,
    eps = eps,
    intercept = intercept
  )
}

# Run subgroup LMTP TMLE ---------------------------------------------------

#' Estimate subgroup-specific LMTP means with a simultaneous TMLE
#'
#' @description
#' Vector-based subgroup TMLE for subgroup-specific LMTP means
#'
#' \deqn{
#'   \psi_j = E[Y(\bar A^d) \mid V \in \mathcal A_j],
#'   \quad j = 1,\ldots,J.
#' }
#'
#' This function mirrors `run_tmle_for_LMTP()`, but targets a vector of subgroup
#' means simultaneously. The plug-in estimate uses targeted regressions
#' `m_star`; the estimated EIF used for variance uses untargeted regressions
#' `m_init`.
#'
#' The subgroup clever covariate at time `t` is
#'
#' \deqn{
#'   H_{t,j}^{obs}(O_i)
#'   =
#'   \frac{1\{i \in \mathcal A_j\}}{\hat P(\mathcal A_j)}
#'   \omega_t(O_i).
#' }
#'
#' The policy-evaluation clever covariate is
#'
#' \deqn{
#'   H_{t,j}^{d}(O_i)
#'   =
#'   \frac{1\{i \in \mathcal A_j\}}{\hat P(\mathcal A_j)}
#'   \omega_{t-1}(O_i) r_t(A_t^d,H_t).
#' }
#'
#' @param ds An `LMTPData` object.
#' @param nuisance_factory An `LMTPNuisanceFactory` object.
#' @param subgroup_funs A function or named list of functions taking `ds$df`
#'   and returning 0/1 subgroup indicators.
#' @param fluctuation An `LMTPFluctuationSubmodel` object.
#' @param alpha Wald interval level.
#'
#' @return An `LMTPFit` with vector-valued estimates and an EIF matrix.
#' @export
run_subgroup_tmle_for_LMTP <- function(ds,
                                       nuisance_factory,
                                       subgroup_funs,
                                       fluctuation = LMTPFluctuationSubmodel$new(),
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

  policy_seq <- nuisance_factory$policy_seq

  if (!inherits(policy_seq, "LMTPPolicySequence")) {
    stop("`nuisance_factory$policy_seq` must inherit from `LMTPPolicySequence`.")
  }

  policy_seq$validate_against_data(ds)

  tau <- ds$tau()
  n <- ds$n

  # Subgroup matrix --------------------------------------------------------

  subgroup_mat <- .make_subgroup_matrix(ds, subgroup_funs)

  if (!is.data.frame(subgroup_mat)) {
    subgroup_mat <- as.data.frame(subgroup_mat, check.names = FALSE)
  }

  if (nrow(subgroup_mat) != n) {
    stop("`subgroup_funs` must return vectors of length `ds$n`.")
  }

  subgroup_mat[] <- lapply(subgroup_mat, as.numeric)

  if (!all(vapply(subgroup_mat, function(x) all(is.finite(x)), logical(1)))) {
    stop("Subgroup indicators must be finite.")
  }

  if (!all(vapply(subgroup_mat, function(x) all(x %in% c(0, 1)), logical(1)))) {
    stop("Each subgroup function must return 0/1 indicators.")
  }

  p_subgroup <- colMeans(subgroup_mat)

  if (any(p_subgroup <= 0)) {
    bad <- colnames(subgroup_mat)[p_subgroup <= 0]
    stop(
      "The following subgroup(s) have zero empirical prevalence: ",
      paste(bad, collapse = ", ")
    )
  }

  subgroup_names <- colnames(subgroup_mat)

  if (is.null(subgroup_names) || any(subgroup_names == "")) {
    subgroup_names <- paste0("subgroup", seq_len(ncol(subgroup_mat)))
    colnames(subgroup_mat) <- subgroup_names
  }

  n_groups <- ncol(subgroup_mat)

  subgroup_weights <- sweep(
    subgroup_mat,
    MARGIN = 2,
    STATS = p_subgroup,
    FUN = "/"
  )

  # 1. Fit ratio nuisance vectors -----------------------------------------

  nuisance_factory$train_ratios(ds)

  # 2. Fit initial, untargeted sequential regressions ----------------------

  m_init_obs <- vector("list", tau)
  m_init_d <- vector("list", tau)

  pseudo_outcome_init <- ds$Y()

  for (t in rev(seq_len(tau))) {
    m_t <- nuisance_factory$train_m_t(
      ds = ds,
      t = t,
      pseudo_outcome = pseudo_outcome_init
    )

    m_init_obs[[t]] <- m_t$m_obs
    m_init_d[[t]] <- m_t$m_d

    if (t > 1L) {
      pseudo_outcome_init <- m_t$m_d
    }
  }

  # 3. Targeted backward recursion ----------------------------------------

  m_star_obs <- vector("list", tau)
  m_star_d <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)
  fluctuation_fits <- vector("list", tau)

  pseudo_outcome_star <- ds$Y()

  for (t in rev(seq_len(tau))) {
    omega_t <- nuisance_factory$omega(t)
    omega_prev_t <- nuisance_factory$omega_prev(t, n = n)
    r_d_t <- nuisance_factory$r_d_preds[[t]]

    H_obs_t <- sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = omega_t,
      FUN = "*"
    )

    H_d_t <- sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = omega_prev_t * r_d_t,
      FUN = "*"
    )

    update_t <- fluctuation$fit_update(
      m_obs = m_init_obs[[t]],
      target = pseudo_outcome_star,
      H_obs = H_obs_t,
      m_d = m_init_d[[t]],
      H_d = H_d_t,
      t = t
    )

    m_star_obs[[t]] <- update_t$m_obs_star
    m_star_d[[t]] <- update_t$m_d_star
    eps[[t]] <- update_t$epsilon
    intercept[[t]] <- update_t$intercept
    fluctuation_fits[[t]] <- update_t$fit

    if (t > 1L) {
      pseudo_outcome_star <- update_t$m_d_star
    }
  }

  # 4. Subgroup plug-in estimates -----------------------------------------

  psi_hat <- colMeans(
    sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = m_star_d[[1L]],
      FUN = "*"
    )
  )

  names(psi_hat) <- subgroup_names

  psi_plugin_init <- colMeans(
    sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = m_init_d[[1L]],
      FUN = "*"
    )
  )

  names(psi_plugin_init) <- subgroup_names

  # 5. Estimated EIF / IC --------------------------------------------------
  #
  # Crucial: use untargeted nuisance regressions m_init, not m_star.
  #
  # Scalar LMTP IC:
  #
  #   ic_i =
  #     sum_t omega_t(i) { next_init_t(i) - m_t(A_t,H_t)(i) }
  #     + m_1(A_1^d,H_1)(i)
  #
  # Subgroup-specific EIF:
  #
  #   D_{j,i} =
  #     [1{i in A_j} / P(A_j)] * ic_i - psi_j.
  #
  # Equivalently:
  #
  #   D_{j,i} =
  #     [1{i in A_j} / P(A_j)]
  #       [
  #         sum_t omega_t(i) { next_init_t(i) - m_t(A_t,H_t)(i) }
  #         + m_1(A_1^d,H_1)(i)
  #       ]
  #     - psi_j.
  #

  ic_scalar_without_centering <- rep(0, n)

  for (t in seq_len(tau)) {
    next_term_untargeted <- if (t == tau) {
      ds$Y()
    } else {
      m_init_d[[t + 1L]]
    }

    ic_scalar_without_centering <- ic_scalar_without_centering +
      nuisance_factory$omega(t) *
      (next_term_untargeted - m_init_obs[[t]])
  }

  ic_scalar_without_centering <- ic_scalar_without_centering + m_init_d[[1L]]

  eif <- sweep(
    subgroup_weights,
    MARGIN = 1,
    STATS = ic_scalar_without_centering,
    FUN = "*"
  )

  eif <- sweep(
    eif,
    MARGIN = 2,
    STATS = psi_hat,
    FUN = "-"
  )

  colnames(eif) <- subgroup_names

  # Optional diagnostic IC using targeted regressions.
  # This is not used for variance by default.
  ic_star_without_centering <- rep(0, n)

  for (t in seq_len(tau)) {
    next_term_targeted <- if (t == tau) {
      ds$Y()
    } else {
      m_star_d[[t + 1L]]
    }

    ic_star_without_centering <- ic_star_without_centering +
      nuisance_factory$omega(t) *
      (next_term_targeted - m_star_obs[[t]])
  }

  ic_star_without_centering <- ic_star_without_centering + m_star_d[[1L]]

  eif_star <- sweep(
    subgroup_weights,
    MARGIN = 1,
    STATS = ic_star_without_centering,
    FUN = "*"
  )

  eif_star <- sweep(
    eif_star,
    MARGIN = 2,
    STATS = psi_hat,
    FUN = "-"
  )

  colnames(eif_star) <- subgroup_names

  # 6. Variance and Wald intervals ----------------------------------------

  var_hat <- apply(eif, 2, stats::var) / n
  se_hat <- sqrt(var_hat)

  ci_hat <- rbind(
    lower = psi_hat - stats::qnorm(1 - alpha / 2) * se_hat,
    upper = psi_hat + stats::qnorm(1 - alpha / 2) * se_hat
  )

  colnames(ci_hat) <- subgroup_names

  # 7. Return fit ----------------------------------------------------------

  LMTPFit$new(
    estimate = psi_hat,
    var = var_hat,
    se = se_hat,
    ci = ci_hat,
    alpha = alpha,
    parameter = "Subgroup E[Y^d]",
    estimator = "Subgroup LMTP TMLE",
    n = n,
    ic = eif,
    eif = eif,

    m_init = list(
      m_obs = m_init_obs,
      m_d = m_init_d
    ),
    m_star = list(
      m_obs = m_star_obs,
      m_d = m_star_d
    ),

    omega = nuisance_factory$omega_preds,
    eps = eps,
    intercept = intercept,
    subgroup_names = subgroup_names
  )
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
