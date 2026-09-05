# Run one-step (AIPW) estimator for LMTP -------------------------------------

#' Estimate an LMTP mean with the one-step (AIPW) estimator
#'
#' @description
#' Estimates the LMTP mean
#'
#' \deqn{
#'   \psi
#'   =
#'   E[Y(\bar A^d)]
#' }
#'
#' with the one-step / AIPW estimator: the plug-in from the initial
#' (untargeted) sequential regressions, plus the empirical mean of the
#' estimated efficient influence function evaluated at those same initial
#' nuisances,
#'
#' \deqn{
#'   \hat\psi_{os}
#'   =
#'   P_n[
#'     \hat\varphi
#'   ],
#'   \qquad
#'   \hat\varphi_i
#'   =
#'   \sum_{t=1}^{\tau}
#'     \omega_t(O_i)
#'     \{
#'       \tilde m_{t+1}(O_i) - m_t(A_t, H_t)(O_i)
#'     \}
#'   +
#'   m_1(A_1^d, H_1)(O_i),
#' }
#'
#' where \eqn{\tilde m_{\tau+1} = Y} and
#' \eqn{\tilde m_{t+1} = m_{t+1}(A_{t+1}^d, H_{t+1})} otherwise.
#'
#' The one-step estimator solves the same efficient-score equation as the
#' TMLE (by construction, \eqn{P_n[\hat D] = 0} at \eqn{\hat\psi_{os}}) but
#' performs no targeting step: the nuisance regressions are never updated.
#' It therefore shares first-order asymptotics with
#' [run_tmle_for_LMTP()] while differing in finite samples, which makes the
#' pair a useful diagnostic: when the two estimators are run on the same data
#' with the same learners, their difference isolates the contribution of the
#' fluctuation/targeting step. In particular, with an identity-link linear
#' fluctuation at \eqn{\tau = 1},
#'
#' \deqn{
#'   \hat\psi_{tmle} - \hat\psi_{os}
#'   =
#'   P_n[\omega (Y - m)]
#'   \left(
#'     \frac{P_n[H^d]}{P_n[\omega^2]} - 1
#'   \right),
#' }
#'
#' so a biased TMLE alongside an unbiased one-step points directly at the
#' heavy-tailed \eqn{P_n[\omega^2]} balance factor rather than at the
#' nuisance estimators.
#'
#' Unlike the TMLE, the one-step plug-in is not guaranteed to respect known
#' outcome bounds.
#'
#' @param ds An `LMTPData` object.
#'
#' @param nuisance_factory An `LMTPNuisanceFactory`.
#'
#' @param alpha Significance level for Wald confidence intervals.
#'
#' @return An `LMTPFit`. The `m_star` slot is `NULL` and `eps` is a list of
#'   zeros, reflecting that no targeting is performed.
#'
#' @seealso [run_tmle_for_LMTP()] for the targeted analog.
#'
#' @export
run_onestep_for_LMTP <- function(ds,
                                 nuisance_factory,
                                 alpha = 0.05) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }
  if (!inherits(nuisance_factory, "LMTPNuisanceFactory")) {
    stop("`nuisance_factory` must inherit from `LMTPNuisanceFactory`.")
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
  # This recursion is identical to the initial step of the TMLE: the
  # one-step and the TMLE share their nuisances exactly when run with the
  # same learners on the same data.
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
  # 3. One-step estimate. ----
  #
  # Construct the uncentered influence representation
  #
  #   phi_i
  #     = sum_t omega_t(O_i)
  #         { next_untargeted_t(O_i) - m_t(A_t,H_t)(O_i) }
  #       + m_1(A_1^d,H_1)(O_i)
  #
  # and estimate
  #
  #   psi_hat = P_n[phi].
  #

  phi <- rep(0, n)

  for (t in seq_len(tau)) {
    next_term_untargeted <- if (t == tau) {
      ds$Y()
    } else {
      m_init_d[[t + 1L]]
    }

    phi <- phi +
      nuisance_factory$omega(t) *
      (next_term_untargeted - m_init_obs[[t]])
  }

  phi <- phi + m_init_d[[1L]]

  psi_hat <- mean(phi)

  #
  # 4. Estimated EIF / IC. ----
  #
  # The one-step solves P_n[D_hat] = 0 by construction:
  #
  #   D_hat(O_i) = phi_i - psi_hat.
  #

  eif <- phi - psi_hat

  var_hat <- stats::var(eif) / n
  se_hat <- sqrt(var_hat)
  ci <- .wald_ci(psi_hat, se_hat, alpha = alpha)

  # 5. Return fitted LMTP object ----

  LMTPFit$new(
    estimate = psi_hat,
    var = var_hat,
    se = se_hat,
    ci = ci,
    alpha = alpha,
    parameter = "E[Y^d]",
    estimator = "One-step (AIPW)",
    n = n,
    ic = eif,
    eif = eif,

    m_init = list(
      m_obs = m_init_obs,
      m_d = m_init_d
    ),
    m_star = NULL,
    omega = nuisance_factory$omega_preds,
    eps = replicate(tau, 0, simplify = FALSE),
    intercept = replicate(tau, 0, simplify = FALSE)
  )
}


# Run subgroup one-step (AIPW) estimator for LMTP -----------------------------

#' Estimate subgroup-specific LMTP means with the one-step (AIPW) estimator
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
#'   j = 1,\ldots,J,
#' }
#'
#' with the one-step / AIPW estimator. Subgroups are specified exactly as in
#' [run_subgroup_tmle_for_LMTP()], via either `subgroup_cols` or
#' `subgroup_funs`; subgroups defined through `subgroup_funs` may overlap.
#'
#' With \eqn{G_{ij} = 1\{i \in \text{subgroup } j\}},
#' \eqn{\hat p_j = P_n G_j}, and \eqn{\hat\varphi} the uncentered scalar LMTP
#' influence representation built from the initial (untargeted) nuisances
#' (see [run_onestep_for_LMTP()]), the estimator is
#'
#' \deqn{
#'   \hat\psi_{os,j}
#'   =
#'   P_n\left[
#'     \frac{G_j}{\hat p_j}
#'     \hat\varphi
#'   \right],
#' }
#'
#' i.e. the empirical mean of \eqn{\hat\varphi} within subgroup j, with
#' subgroup EIF
#'
#' \deqn{
#'   \hat D_j(O_i)
#'   =
#'   \frac{G_{ij}}{\hat p_j}
#'   \{
#'     \hat\varphi_i - \hat\psi_{os,j}
#'   \}.
#' }
#'
#' Each subgroup score \eqn{P_n[\hat D_j] = 0} holds by construction, with no
#' targeting step. Run alongside [run_subgroup_tmle_for_LMTP()] with the same
#' learners on the same data, the difference between the two fits isolates
#' the contribution of the simultaneous vector fluctuation, which is the
#' recommended first diagnostic when subgroup TMLE estimates look biased.
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
#' @param alpha Significance level for Wald confidence intervals.
#'
#' @return An `LMTPFit` with vector-valued estimates and an n x J EIF matrix.
#'   The `m_star` slot is `NULL` and `eps` is a list of zero vectors,
#'   reflecting that no targeting is performed.
#'
#' @seealso [run_subgroup_tmle_for_LMTP()] for the targeted analog.
#'
#' @export
run_subgroup_onestep_for_LMTP <- function(
    ds,
    nuisance_factory,
    subgroup_cols = NULL,
    subgroup_funs = NULL,
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
  # Identical construction and validation to run_subgroup_tmle_for_LMTP():
  # an n x J matrix G with G[i, j] = 1 if observation i belongs to
  # subgroup j.

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
  #   subgroup_weights[i, j] = G[i, j] / p_hat_j.

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
  # Identical recursion to the TMLE runners:
  #
  #   pseudo_tau = Y
  #   pseudo_{t-1} = m_t(A_t^d, H_t).

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


  # 3. Construct the uncentered scalar influence representation ----
  #
  #   phi_i
  #     = sum_t omega_t(i) { next_init_t(i) - m_t(A_t, H_t)(i) }
  #       + m_1(A_1^d, H_1)(i).

  phi <- rep(0, n)

  for (t in seq_len(tau)) {

    next_term_untargeted <- if (t == tau) {
      ds$Y()
    } else {
      m_init_d[[t + 1L]]
    }

    phi <- phi +
      nuisance_factory$omega(t) *
      (
        next_term_untargeted -
          m_init_obs[[t]]
      )
  }

  phi <- phi + m_init_d[[1L]]


  # 4. Compute subgroup one-step estimates ----
  #
  # For subgroup j,
  #
  #   psi_hat_j = P_n[ G_j / p_hat_j * phi ],
  #
  # the empirical mean of phi among observations in subgroup j.

  psi_hat <- colMeans(
    sweep(
      subgroup_weights,
      MARGIN = 1,
      STATS = phi,
      FUN = "*"
    )
  )

  names(psi_hat) <- subgroup_names


  # 5. Compute subgroup efficient influence functions ----
  #
  #   D_j(O_i) = G_ij / p_hat_j * { phi_i - psi_hat_j }.
  #
  # The centering by psi_hat_j occurs INSIDE the subgroup weight, exactly as
  # in run_subgroup_tmle_for_LMTP(). P_n[D_j] = 0 for each j by construction.

  centered_ic <- matrix(
    phi,
    nrow = n,
    ncol = n_groups
  )

  centered_ic <- sweep(
    centered_ic,
    MARGIN = 2,
    STATS = psi_hat,
    FUN = "-"
  )

  eif <- as.matrix(
    subgroup_weights
  ) * centered_ic

  colnames(eif) <- subgroup_names


  # 6. Compute variances and pointwise Wald intervals ----

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
    estimator = "Subgroup LMTP one-step (AIPW)",
    n = n,
    ic = eif,
    eif = eif,

    m_init = list(
      m_obs = m_init_obs,
      m_d = m_init_d
    ),

    m_star = NULL,
    omega = nuisance_factory$omega_preds,
    eps = replicate(
      tau,
      stats::setNames(rep(0, n_groups), subgroup_names),
      simplify = FALSE
    ),
    intercept = replicate(tau, 0, simplify = FALSE),
    subgroup_names = subgroup_names
  )
}
