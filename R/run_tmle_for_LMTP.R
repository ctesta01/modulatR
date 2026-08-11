
# this file handles the core runners for estimating E[Y^d]  and E[Y^d | V = v]

#' Estimate a scalar TMLE for an LMTP mean
#'
#' @description
#' TMLE runner for estimating
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
  # 1. Fit ratio nuisance vectors. ----
  #
  # Stores:
  #   r_obs[[t]]     = r_t(A_t, H_t)
  #   r_d[[t]]       = r_t(A_t^d, H_t)
  #   omega_obs[[t]] = prod_{s <= t} r_s(A_s, H_s)
  #

  nuisance_factory$train_ratios(ds)

  #
  # 2. Fit initial, untargeted sequential regressions. ----
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
  # 3. Targeted backward recursion. ----
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
  # 4. TMLE plug-in estimate. ----
  #
  # Uses targeted m_1 evaluated at the modified treatment value.
  #

  psi_hat <- mean(m_star_d[[1L]])

  #
  # 5. Estimated EIF / IC. ----
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

  # 6. Return fitted LMTP object ----

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
#' Estimates the vector of subgroup-specific LMTP means
#'
#' \deqn{
#'   \psi_j
#'   =
#'   E[
#'     Y(\bar A^d)
#'     \mid
#'     V = v_j
#'   ],
#'   \qquad
#'   j = 1,\ldots,J.
#' }
#'
#' Subgroups can be specified in one of two ways:
#'
#' \itemize{
#'   \item
#'   `subgroup_cols`: One or more categorical
#'   baseline columns are supplied, and each observed combination of their
#'   values defines one subgroup.
#'
#'   \item
#'   `subgroup_funs`: Or, supply a function or named list
#'   of functions mapping `ds$df` to 0/1 subgroup indicators. This allows
#'   arbitrary and potentially overlapping subgroups.
#' }
#'
#' Exactly one of `subgroup_cols` and `subgroup_funs` must be supplied.
#'
#' Let
#'
#' \deqn{
#'   G_{ij}
#'   =
#'   1\{i \text{ belongs to subgroup } j\},
#' }
#'
#' and let
#'
#' \deqn{
#'   \hat p_j
#'   =
#'   P_n G_j.
#' }
#'
#' The corresponding normalized subgroup weight is
#'
#' \deqn{
#'   W_{ij}
#'   =
#'   \frac{G_{ij}}{\hat p_j}.
#' }
#'
#' These weights are used in both the simultaneous targeting step and the
#' subgroup efficient influence functions.
#'
#' @param ds An `LMTPData` object.
#'
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#'
#' @param subgroup_cols Optional character vector naming categorical baseline
#'   columns in `ds$W_cols`. Each observed combination of these columns defines
#'   one subgroup.
#'
#' @param subgroup_funs Optional function or named list of functions taking
#'   `ds$df` and returning 0/1 subgroup indicators. Intended as an advanced
#'   interface for arbitrary subgroup definitions. Subgroups may overlap.
#'
#' @param fluctuation An `LMTPFluctuationSubmodel`.
#'
#' @param alpha Significance level for Wald confidence intervals.
#'
#' @return An `LMTPFit` with vector-valued estimates and an n x J EIF matrix.
#'
#' @export
run_subgroup_tmle_for_LMTP <- function(
    ds,
    nuisance_factory,
    subgroup_cols = NULL,
    subgroup_funs = NULL,
    fluctuation = LMTPFluctuationSubmodel$new(),
    alpha = 0.05
) {

  # Validate inputs ----

  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) {
    stop(
      "`nuisance_factory` must inherit from ",
      "`LMTPNuisanceFactory`."
    )
  }

  if (!inherits(fluctuation, "LMTPFluctuationSubmodel")) {
    stop(
      "`fluctuation` must inherit from ",
      "`LMTPFluctuationSubmodel`."
    )
  }

  has_cols <- !is.null(subgroup_cols)
  has_funs <- !is.null(subgroup_funs)

  if (has_cols == has_funs) {
    stop(
      "Supply exactly one of `subgroup_cols` or `subgroup_funs`."
    )
  }

  policy_seq <- nuisance_factory$policy_seq

  if (!inherits(policy_seq, "LMTPPolicySequence")) {
    stop(
      "`nuisance_factory$policy_seq` must inherit from ",
      "`LMTPPolicySequence`."
    )
  }

  policy_seq$validate_against_data(ds)

  tau <- ds$tau()
  n <- ds$n


  # Construct subgroup indicator matrix ----
  #
  # Represent the J subgroups using an n x J matrix G:
  #
  #   G[i, j] = 1 if observation i belongs to subgroup j,
  #             0 otherwise.
  #
  # When `subgroup_cols` is used, each distinct observed combination of the
  # conditioning variables defines one mutually exclusive subgroup.
  #
  # When `subgroup_funs` is used, each user-supplied function defines an
  # indicator column, so groups may overlap.

  subgroup_mat <- if (has_cols) {

    .make_subgroup_indicators(
      ds = ds,
      subgroup_cols = subgroup_cols
    )

  } else {

    .make_subgroup_matrix(
      ds = ds,
      subgroup_funs = subgroup_funs
    )
  }

  if (!is.data.frame(subgroup_mat)) {
    subgroup_mat <- as.data.frame(
      subgroup_mat,
      check.names = FALSE
    )
  }

  if (nrow(subgroup_mat) != n) {
    stop(
      "The subgroup specification must produce ",
      "one indicator value per observation."
    )
  }

  subgroup_mat[] <- lapply(
    subgroup_mat,
    as.numeric
  )

  if (!all(
    vapply(subgroup_mat, function(x)
      all(is.finite(x)), logical(1))
  )) {
    stop("Subgroup indicators must be finite.")
  }

  if (!all(
    vapply(subgroup_mat, function(x)
      all(x %in% c(0, 1)), logical(1))
  )) {
    stop(
      "Each subgroup indicator must contain only 0 and 1."
    )
  }

  subgroup_names <- colnames(subgroup_mat)

  if (
    is.null(subgroup_names) ||
    any(subgroup_names == "")
  ) {
    subgroup_names <- paste0(
      "subgroup",
      seq_len(ncol(subgroup_mat))
    )

    colnames(subgroup_mat) <- subgroup_names
  }

  n_groups <- ncol(subgroup_mat)


  # Compute normalized subgroup weights ----
  #
  # Estimate the subgroup prevalence
  #
  #   p_hat_j = P_n(G_j = 1).
  #
  # We repeatedly use the identity
  #
  #   E[X | G_j = 1]
  #   =
  #   E[
  #     G_j / P(G_j = 1) * X
  #   ].
  #
  # Therefore define
  #
  #   subgroup_weights[i, j]
  #   =
  #   G[i, j] / p_hat_j.
  #
  # Each column is zero outside subgroup j and equal to 1 / p_hat_j inside it.

  p_subgroup <- colMeans(
    subgroup_mat
  )

  if (any(p_subgroup <= 0)) {
    bad <- subgroup_names[
      p_subgroup <= 0
    ]

    stop(
      "The following subgroup(s) have zero empirical prevalence: ",
      paste(bad, collapse = ", ")
    )
  }

  subgroup_weights <- sweep(
    subgroup_mat,
    MARGIN = 2,
    STATS = p_subgroup,
    FUN = "/"
  )


  # 1. Fit treatment-density ratios ----
  nuisance_factory$train_ratios(ds)

  # 2. Fit initial sequential regressions ----
  #
  # Work backwards from the terminal outcome:
  #
  #   pseudo_tau = Y
  #
  # and, recursively,
  #
  #   pseudo_{t-1}
  #   =
  #   m_t(A_t^d, H_t).
  #
  # Store each initial regression at both the observed and shifted treatment:
  #
  #   m_init_obs[[t]] = m_t(A_t, H_t)
  #   m_init_d[[t]]   = m_t(A_t^d, H_t)

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


  # 3. Perform simultaneous subgroup targeting ----
  #
  # At time t and for subgroup j, the observed-data clever covariate is
  #
  #   H_obs[i, j]
  #   =
  #   G[i, j] / p_hat_j
  #   * omega_t(i).
  #
  # The clever covariate used to evaluate the updated regression under the
  # modified treatment policy is
  #
  #   H_d[i, j]
  #   =
  #   G[i, j] / p_hat_j
  #   * omega_{t-1}(i)
  #   * r_t(A_t^d, H_t).
  #
  # The `sweep(..., MARGIN = 1, ...)` operations below multiply each row i
  # of the subgroup-weight matrix by the corresponding scalar longitudinal
  # weight for observation i.

  m_star_obs <- vector("list", tau)
  m_star_d <- vector("list", tau)
  eps <- vector("list", tau)
  intercept <- vector("list", tau)
  pseudo_outcome_star <- ds$Y()

  for (t in rev(seq_len(tau))) {

    omega_t <- nuisance_factory$omega(t)

    omega_prev_t <- nuisance_factory$omega_prev(t, n = n)

    r_d_t <- nuisance_factory$r_d_preds[[t]]

    # H_obs[i, j] = G[i, j] / p_hat_j * omega_t(i)
    H_obs_t <- sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = omega_t,
      FUN = "*"
    )

    # H_d[i, j] = G[i, j] / p_hat_j  * omega_{t-1}(i)  * r_t(A_t^d, H_t)
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

    if (t > 1L) {
      pseudo_outcome_star <- update_t$m_d_star
    }
  }


  # 4. Compute subgroup plug-in estimates ----
  #
  # For subgroup j,
  #
  #   psi_hat_j
  #   =
  #   P_n[
  #     G_j / p_hat_j
  #     * m_1^d,star
  #   ].
  #
  # Since p_hat_j = n_j / n, this is exactly the empirical mean of
  # m_1^d,star among observations belonging to subgroup j.

  psi_hat <- colMeans(
    sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = m_star_d[[1L]],
      FUN = "*"
    )
  )

  names(psi_hat) <- subgroup_names


  # 5. Compute subgroup efficient influence functions ----
  #
  # First construct the uncentered scalar LMTP influence representation
  #
  #   phi_i
  #   =
  #   sum_t
  #     omega_t(i)
  #     {
  #       next_init_t(i)
  #       -
  #       m_t(A_t, H_t)(i)
  #     }
  #   +
  #   m_1(A_1^d, H_1)(i).
  #
  # The subgroup-specific EIF is then
  #
  #   D_j(O_i)
  #   =
  #   G_ij / p_j
  #   {
  #     phi_i - psi_j
  #   }.
  #
  # The centering by psi_j occurs INSIDE the subgroup weight.
  #
  # It is not
  #
  #   G_ij / p_j * phi_i - psi_j.

  ic_scalar_without_centering <- rep(
    0,
    n
  )

  for (t in seq_len(tau)) {

    next_term_untargeted <- if (t == tau) {
      ds$Y()
    } else {
      m_init_d[[t + 1L]]
    }

    ic_scalar_without_centering <-
      ic_scalar_without_centering +
      nuisance_factory$omega(t) *
      (
        next_term_untargeted -
          m_init_obs[[t]]
      )
  }

  ic_scalar_without_centering <-
    ic_scalar_without_centering +
    m_init_d[[1L]]

  # Construct an n x J matrix whose (i, j) entry is phi_i - psi_hat_j.
  centered_ic <- matrix(
    ic_scalar_without_centering,
    nrow = n,
    ncol = n_groups
  )

  centered_ic <- sweep(
    centered_ic,
    MARGIN = 2,
    STATS = psi_hat,
    FUN = "-"
  )

  # EIF[i, j]  =  G[i, j] / p_hat_j * {phi_i - psi_hat_j}
  eif <- as.matrix(
    subgroup_weights
  ) * centered_ic

  colnames(eif) <- subgroup_names


  # 6. Compute variances and pointwise Wald intervals ----
  #
  # For subgroup j,
  #
  #   Var(psi_hat_j)
  #   approximately equals
  #   Var(D_j) / n.
  #
  # These will be pointwise confidence intervals. Simultaneous inference can
  # subsequently be added using the multiplier bootstrap.

  var_hat <- apply(
    eif,
    MARGIN = 2,
    FUN = stats::var
  ) / n

  se_hat <- sqrt(var_hat)

  z <- stats::qnorm(1 - alpha / 2)

  ci_hat <- rbind(
    lower = psi_hat - z * se_hat,
    upper = psi_hat + z * se_hat
  )

  colnames(ci_hat) <- subgroup_names


  # 7. Return fitted subgroup LMTP object ----

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
