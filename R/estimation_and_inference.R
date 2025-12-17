# TMLE / G-computation / IPW for LMTP -----------------------------------------

#' LMTP estimation via TMLE, outcome regression, or IPW
#'
#' @description
#' Core workhorse to estimate the marginal LMTP parameter
#' `psi = E[Y^d]` under a longitudinal modified treatment policy. This
#' function performs a sequential outcome regression, optionally applies
#' a TMLE fluctuation using the LMTP clever covariate, and computes the
#' efficient influence curve using the *unfluctuated* outcome regression
#' trajectory `Q0_t` (as recommended for variance estimation in
#' longitudinal TMLE; see discussion in Díaz et al. 2023 and related
#' work). The plug-in estimate is still based on the fluctuated
#' trajectory when `method = "tmle"`.
#'
#' @param ds `LMTP_Data_Struct` instance.
#' @param policy_seq `LMTPPolicySequence` describing the LMTP.
#' @param learners_Q List (or single spec) of learner specifications for
#'   sequential outcome regressions.
#' @param learners_g_factory An `LMTPNuisanceFactory` instance (already
#'   configured, but not yet trained).
#' @param fml_Q Optional formula or list of formulas for the outcome
#'   regressions.
#' @param outcome_link Either `"identity"` or `"logit"`.
#' @param bounds Numeric vector of length 2 giving outcome bounds when
#'   `outcome_link = "logit"`.
#' @param maxit Maximum number of outer TMLE iterations.
#' @param eps_tol Convergence tolerance for the fluctuation (ignored when
#'   `maxit = 1`).
#' @param repeat_lnrs Logical; if `TRUE`, recycle `learners_Q` across
#'   time points when only one is supplied.
#' @param method One of `"tmle"`, `"outcome"` (sequential G-computation
#'   without targeting), or `"ipw"`.
#'
#' @return An object of class `tmle_lmtp` with components:
#'   * `psi`: point estimate of `E[Y^d]`.
#'   * `se`: standard error.
#'   * `ci95`: 95\% Wald confidence interval.
#'   * `ic`: estimated influence curve values.
#'   * `Q_star`: list of fluctuated regression functions `Q_t^*`.
#'   * `Q0`: list of unfluctuated regression functions `Q0_t`.
#'   * `Kprov`: the K-provider used.
#'
#' @export
fit_tmle_for_LMTP <- function(
    ds,
    policy_seq,
    learners_Q,
    learners_g_factory,
    fml_Q = NULL,
    outcome_link = c("identity", "logit"),
    bounds = c(0, 1),
    maxit = 1,
    eps_tol = 1e-6,
    repeat_lnrs = TRUE,
    method = c("tmle", "outcome", "ipw")
) {
  method <- match.arg(method)
  outcome_link <- match.arg(outcome_link)
  tau <- ds$tau()

  # 1) Fit g and r once
  learners_g_factory$train(ds)

  # 2) Helpers from the factory
  Kprov <- learners_g_factory$k_provider(ds)
  learners_Q <- .rep_if_needed(learners_Q, ds$tau(), repeat_lnrs)

  Q_trainer_fun <- learners_g_factory$q_trainer(
    learners_Q = if (length(learners_Q) == 1) rep(list(learners_Q[[1]]), tau) else learners_Q,
    fml_Q = fml_Q,
    repeat_fmls = is.null(fml_Q) || length(fml_Q) == 1
  )

  fluct <- switch(outcome_link,
                  identity = LMTPFluctuationIdentity$new(),
                  logit    = LMTPFluctuationLogit$new(bounds = bounds))

  # 3) Initialize storage for Q_star and Q0 closures
  Q_star <- vector("list", tau)
  Q0_list <- vector("list", tau)

  # 4) Outer iterate
  iter <- 1; eps <- 1e6
  while (iter <= maxit && eps > eps_tol) {

    tilde_next <- ds$Y()

    for (t in tau:1) {
      Ht <- ds$H(t)

      # Train initial Q0_t on pseudo-outcome
      Q0_t <- Q_trainer_fun(ds, t, tilde_next)
      Q0_list[[t]] <- Q0_t

      # target for epsilon fit
      target_vec <- if (t < tau) {
        Hnext <- ds$H(t + 1)
        Astar_next <- policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext)
        if (!is.null(Q_star[[t + 1]])) {
          Q_star[[t + 1]](Astar_next, Hnext)
        } else {
          Q0_list[[t + 1]](Astar_next, Hnext)
        }
      } else {
        ds$Y()
      }

      # clever covariate and offset evals on observed data
      K_obs_t <- Kprov$K_obs(t)
      Qvec <- Q0_t(ds$A(t), Ht)

      if (method == "tmle") {
        up <- fluct$fit_update(
          Q0_fun_t = Q0_t,
          Qvec = Qvec,
          target_vec = target_vec,
          K_obs_t = K_obs_t,
          Kprov = Kprov,
          t = t)

        eps <- up$epsilon
        Q_star[[t]] <- up$wrap
      } else {
        Q_star[[t]] <- Q0_t
        eps <- 0
      }

      # carry upstream the pseudo-outcome under the policy
      Astar_t <- policy_seq$apply_policy_t(t, ds$A(t), Ht)
      tilde_next <- Q_star[[t]](Astar_t, Ht)
    }

    iter <- iter + 1
  }

  # 5) Estimate psi
  H1 <- ds$H(1)
  Astar1 <- policy_seq$apply_policy_t(1, ds$A(1), H1)
  if (method %in% c("tmle", "outcome")) {
    psi_hat <- mean(Q_star[[1]](Astar1, H1))
  } else if (method == "ipw") {
    psi_hat <- mean(Kprov$K_obs(tau) * ds$Y())
  }

  # 6) Influence curve using *unfluctuated* Q0_t
  if (method == "tmle") {
    n <- nrow(ds$data)
    phi <- rep(0, n)
    for (t in seq_len(tau)) {
      Ht <- ds$H(t)
      mt_obs <- Q0_list[[t]](ds$A(t), Ht)
      mnext <- if (t < tau) {
        Hnext <- ds$H(t + 1)
        Astar_next <- policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext)
        Q0_list[[t + 1]](Astar_next, Hnext)
      } else ds$Y()
      w <- Kprov$K_obs(t)
      phi <- phi + w * (mnext - mt_obs)
    }
    plug <- Q_star[[1]](Astar1, H1)
    ic <- phi + plug - psi_hat
    se <- sqrt(stats::var(ic) / n)
    ci95 <- c(psi_hat - 1.96 * se, psi_hat + 1.96 * se)
  } else {
    se <- NA; ci95 <- c(NA, NA); ic <- NA
  }

  structure(list(
    psi = psi_hat,
    se = se,
    ci95 = ci95,
    ic = ic,
    Q_star = Q_star,
    Q0 = Q0_list,
    Kprov = Kprov
  ), class = c("tmle_lmtp", "list"))
}

# Heterogeneous LMTP effects ---------------------------------------------------

#' Heterogeneous LMTP effects for baseline subgroups
#'
#' @description
#' Estimate subgroup-specific LMTP effects
#'  \eqn{\theta_v = E[Y^d | V = v]} for a set of baseline subgroups defined by
#' baseline covariates `V ⊂ W`. The procedure implements the joint TMLE
#' described in the prompt: a single set of nuisance fits and a
#' multi-dimensional clever covariate that targets all subgroup effects
#' simultaneously, following the strategy of Wei et al. (2023) adapted to
#' the LMTP setting.
#'
#' Subgroups are defined as the unique combinations of columns in `V`.
#' For each subgroup level `v`, the target parameter is
#' `theta_v = E[Y^d | V = v]` and the EIF is
#'
#' \deqn{D_v^*(O) = \frac{I(V=v)}{p_v}\left\{\sum_{t=1}^\tau
#'   \omega_t( m_{t+1}(A_t^d, H_t) - E[m_{t+1}(A_{t+1}, H_{t+1}) | A_t, H_t]) + m_1(A_1^d, L_1) - \theta_v \right\},}
#'
#' with \eqn{p_v = P(V = v)} and \eqn{\omega_t = \prod_{s=1}^t r_s(A_s, H_s)}.
#' In this implementation, the TMLE produces a targeted (fluctuated) substitution estimator \eqn{m_1^*},
#' and uses an unfluctuated substitution estimators \eqn{m_1} in the influence curve used for variance
#' calculations and inference.
#'
#' @param ds `LMTP_Data_Struct` instance.
#' @param policy_seq `LMTPPolicySequence` describing the LMTP.
#' @param learners_Q List (or single spec) of learner specifications for
#'   the outcome regressions.
#' @param learners_g_factory An `LMTPNuisanceFactory` instance (already
#'   configured, but not yet trained).
#' @param V Character vector of baseline column names in `ds$W()` used to
#'   define subgroups. All unique combinations of these columns define
#'   the set \eqn{\mathcal V} over which subgroup effects are estimated.
#' @param fml_Q Optional formula or list of formulas for the outcome
#'   regressions.
#' @param outcome_link Either `"identity"` or `"logit"`.
#' @param bounds Bounds for the outcome when `outcome_link = "logit"`.
#' @param maxit Currently ignored (one-step TMLE); included for future
#'   compatibility.
#' @param eps_tol Currently ignored; included for compatibility.
#' @param repeat_lnrs Logical; if `TRUE`, recycle `learners_Q` across
#'   time points when only one is supplied.
#' @param method One of `"tmle"`, `"outcome"`, or `"ipw"`. When
#'   `method = "ipw"`, only inverse probability weighted estimators are
#'   computed, without outcome regression.
#' @param stabilization Character string controlling potential future
#'   stabilization of the subgroup clever covariates / weights. Currently
#'   only `"naive"` is implemented; `"dual"` is reserved for a future
#'   implementation of the Lagrangian-dual approach of Wei et al. (2023).
#'
#' @return An object of class `tmle_lmtp_heterogeneous` with components:
#'   * `theta`: vector of subgroup-specific estimates.
#'   * `se`: standard errors for each subgroup.
#'   * `ci95`: matrix of 95\% Wald intervals (columns `lower`, `upper`).
#'   * `ic`: `n x J` matrix of influence curve values, with columns
#'      corresponding to subgroups.
#'   * `V_levels`: factor levels of the subgroup variable (unique
#'      combinations of `V`).
#'   * `p_v`: empirical subgroup proportions.
#'   * `method`: estimation method used.
#'   * `subgroup_vars`: names of the baseline variables that define the
#'     subgroups.
#'
#' @export
fit_lmtp_heterogeneous <- function(
    ds,
    policy_seq,
    learners_Q,
    learners_g_factory,
    V,
    fml_Q = NULL,
    outcome_link = c("identity", "logit"),
    bounds = c(0, 1),
    maxit = 1,
    eps_tol = 1e-6,
    repeat_lnrs = TRUE,
    method = c("tmle", "outcome", "ipw"),
    stabilization = c("naive", "dual")) {

  outcome_link <- match.arg(outcome_link)
  method <- match.arg(method)
  stabilization <- match.arg(stabilization)

  if (stabilization != "naive") {
    stop("stabilization = 'dual' is not yet implemented; only 'naive' is currently supported.")
  }

  tau <- ds$tau()
  n <- nrow(ds$data)

  # Subgroup variable V: factor of unique baseline combinations
  W_df <- ds$W()
  if (is.null(W_df)) stop("No baseline W_cols supplied in LMTP_Data_Struct; cannot define V.")
  if (!all(V %in% colnames(W_df))) {
    missing_V <- setdiff(V, colnames(W_df))
    stop("Columns not found in baseline W: ", paste(missing_V, collapse = ", "))
  }
  V_df <- W_df[, V, drop = FALSE]
  V_str <- interaction(V_df, drop = TRUE, lex.order = TRUE)
  V_levels <- levels(V_str)
  J <- length(V_levels)

  p_v <- prop.table(table(V_str))
  p_v_vec <- as.numeric(p_v)[match(V_levels, names(p_v))]

  # Train nuisance functions
  learners_g_factory$train(ds)
  Kprov <- learners_g_factory$k_provider(ds)
  learners_Q <- .rep_if_needed(learners_Q, ds$tau(), repeat_lnrs)

  Q_trainer_fun <- learners_g_factory$q_trainer(
    learners_Q = if (length(learners_Q) == 1) rep(list(learners_Q[[1]]), tau) else learners_Q,
    fml_Q = fml_Q,
    repeat_fmls = is.null(fml_Q) || length(fml_Q) == 1
  )

  # Precompute omega_t and subgroup clever covariates K_{t,v}
  omega_list <- lapply(seq_len(tau), function(t) Kprov$K_obs(t))
  K_list <- vector("list", tau)
  for (t in seq_len(tau)) {
    omega_t <- omega_list[[t]]
    K_t <- matrix(0, nrow = n, ncol = J)
    colnames(K_t) <- paste0("K", seq_len(J))
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      if (length(idx) > 0) {
        K_t[idx, j] <- omega_t[idx] / p_v_vec[j]
      }
    }
    K_list[[t]] <- K_t
  }

  # Storage for Q0 and Q* trajectories (vectors evaluated on the sample)
  Q0_obs_list <- vector("list", tau)
  Q0_d_list   <- vector("list", tau)
  Qstar_obs_list <- vector("list", tau)
  Qstar_d_list   <- vector("list", tau)

  # One-step TMLE / G-computation backward recursion
  tilde_next <- ds$Y()

  for (t in tau:1) {
    Ht <- ds$H(t)
    At <- ds$A(t)
    Astar_t <- policy_seq$apply_policy_t(t, At, Ht)

    Q0_t <- Q_trainer_fun(ds, t, tilde_next)
    Q0_obs <- Q0_t(At, Ht)
    Q0_d   <- Q0_t(Astar_t, Ht)

    Q0_obs_list[[t]] <- Q0_obs
    Q0_d_list[[t]]   <- Q0_d

    if (method == "ipw") {
      # no targeting / outcome regression used only for diagnostics
      Qstar_obs <- Q0_obs
      Qstar_d   <- Q0_d
    } else {
      K_t <- K_list[[t]]
      target_vec <- tilde_next

      if (outcome_link == "identity") {
        df <- data.frame(target = target_vec,
                         offset = Q0_obs,
                         K_t)
        fit <- stats::glm(target ~ -1 + offset(offset) + ., data = df, family = gaussian())
        eps <- coef(fit)
        eps[is.na(eps)] <- 0
        eps <- eps[setdiff(names(eps), "offset")]  # keep only K coefficients
        eps_vec <- as.numeric(eps)
        eps_vec <- if (length(eps_vec) < J) c(eps_vec, rep(0, J - length(eps_vec))) else eps_vec
        Qstar_obs <- as.numeric(Q0_obs + K_t %*% eps_vec)
        Qstar_d   <- as.numeric(Q0_d   + K_t %*% eps_vec)
      } else {
        # logit link with bounds
        a <- bounds[1]; b <- bounds[2]
        to01 <- function(x) (x - a) / (b - a)
        from01 <- function(z) a + (b - a) * z
        b01 <- function(z) pmin(pmax(z, 1e-6), 1 - 1e-6)

        Y01 <- to01(target_vec)
        off <- qlogis(b01(to01(Q0_obs)))
        df <- data.frame(Y01 = Y01,
                         offset = off,
                         K_t)
        fit <- stats::glm(Y01 ~ -1 + offset(offset) + ., data = df, family = binomial())
        eps <- coef(fit)
        eps[is.na(eps)] <- 0
        eps <- eps[setdiff(names(eps), "offset")]  # only K coefficients
        eps_vec <- as.numeric(eps)
        eps_vec <- if (length(eps_vec) < J) c(eps_vec, rep(0, J - length(eps_vec))) else eps_vec

        z_obs <- off + as.numeric(K_t %*% eps_vec)
        z_d   <- qlogis(b01(to01(Q0_d))) + as.numeric(K_t %*% eps_vec)
        Qstar_obs <- from01(plogis(z_obs))
        Qstar_d   <- from01(plogis(z_d))
      }
    }

    Qstar_obs_list[[t]] <- Qstar_obs
    Qstar_d_list[[t]]   <- Qstar_d

    # Update pseudo-outcome for next step
    tilde_next <- if (method == "tmle") Qstar_d else Q0_d
  }

  # Plug-in estimators for each subgroup
  theta <- numeric(J)
  names(theta) <- V_levels

  if (method %in% c("tmle", "outcome")) {
    Q1_d <- if (method == "tmle") Qstar_d_list[[1]] else Q0_d_list[[1]]
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      theta[j] <- mean(Q1_d[idx] / p_v_vec[j])
    }
  } else if (method == "ipw") {
    omega_tau <- omega_list[[tau]]
    Y <- ds$Y()
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      theta[j] <- mean(omega_tau[idx] * Y[idx] / p_v_vec[j])
    }
  }

  # Influence curve for tmle / outcome using unfluctuated Q0_t
  if (method %in% c("tmle", "outcome")) {
    phi <- rep(0, n)
    for (t in seq_len(tau)) {
      mt_obs <- Q0_obs_list[[t]]
      mnext <- if (t < tau) Q0_d_list[[t + 1]] else ds$Y()
      phi <- phi + omega_list[[t]] * (mnext - mt_obs)
    }
    plug <- Q0_d_list[[1]]
    ic_mat <- matrix(0, nrow = n, ncol = J)
    colnames(ic_mat) <- V_levels
    for (j in seq_len(J)) {
      idx <- which(V_str == V_levels[j])
      ic_v <- phi + plug - theta[j]
      ic_mat[, j] <- 0
      ic_mat[idx, j] <- ic_v[idx] / p_v_vec[j]
    }
    se <- sqrt(colSums((ic_mat - matrix(colMeans(ic_mat), nrow = n, ncol = J, byrow = TRUE))^2) / (n - 1) / n)
    ci95 <- cbind(lower = theta - 1.96 * se,
                  upper = theta + 1.96 * se)
  } else {
    ic_mat <- matrix(NA_real_, nrow = n, ncol = J)
    colnames(ic_mat) <- V_levels
    se <- rep(NA_real_, J)
    ci95 <- cbind(lower = rep(NA_real_, J), upper = rep(NA_real_, J))
  }

  res <- list(
    theta = theta,
    se = se,
    ci95 = ci95,
    ic = ic_mat,
    V_levels = V_levels,
    p_v = p_v_vec,
    method = method,
    subgroup_vars = V
  )
  class(res) <- c("tmle_lmtp_heterogeneous", "list")
  res
}
