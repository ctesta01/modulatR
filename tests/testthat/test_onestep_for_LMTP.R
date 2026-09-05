# ============================================================================
# Tests for run_onestep_for_LMTP() and run_subgroup_onestep_for_LMTP()
# ============================================================================

library(modulatR)

set.seed(20260831)
`%~%` <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8))
check <- function(cond, msg) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) stop("Test failed: ", msg)
}

gen <- function(n, beta1 = 7.5, pV1 = 0.33) {
  L <- rnorm(n)
  V <- rbinom(n, 1, pV1)
  A <- rnorm(n, mean = L, sd = 1)
  Y <- L + 5 * A + (beta1 - 5) * V * A + rnorm(n)
  data.frame(L = L, V = V, A = A, Y = Y)
}

make_ds <- function(df) {
  modulatR::LMTPData$new(
    data = df, A_cols = "A", L_cols = "L", W_cols = "V", Y_col = "Y"
  )
}

make_factory <- function(ds, m_formula = as.formula("..pseudo_outcome ~ 1")) {
  modulatR::LMTPNuisanceFactory$new(
    policy_seq = modulatR::repeat_policy_over_time(
      mtp = modulatR::mtp_additive_shift(delta = 1, name = "shift1"),
      tau = ds$tau(), name = "shift1_seq"
    ),
    m_learners = modulatR::make_glm_m_learners(
      tau = ds$tau(), family = gaussian(), formula = m_formula
    ),
    g_learners = modulatR::make_nadir_density_g_learners(
      tau = ds$tau(), learner = nadir::lnr_glm_density,
      formula = as.formula("A ~ L")
    )
  )
}

# ---------------------------------------------------------------------------
# 1. Scalar one-step: internal identities
# ---------------------------------------------------------------------------
cat("\n--- 1. run_onestep_for_LMTP ---\n")

df  <- gen(4000)
df1 <- df[df$V == 1, , drop = FALSE]
ds1 <- make_ds(df1)

os <- modulatR::run_onestep_for_LMTP(ds1, make_factory(ds1))

# 1a. Closed form: psi = Pn[m_d] + Pn[omega (Y - m_obs)] at tau = 1
psi_manual <- mean(os$m_init$m_d[[1]]) +
  mean(os$omega[[1]] * (df1$Y - os$m_init$m_obs[[1]]))
check(os$estimate %~% psi_manual,
      "one-step equals Pn[m_d] + Pn[omega (Y - m)]")

# 1b. Score solved by construction
check(mean(os$eif) %~% 0, "Pn[EIF] = 0 at the one-step estimate")

# 1c. SE from the same EIF formula as the TMLE
check(os$se %~% sqrt(stats::var(os$eif) / nrow(df1)),
      "SE matches var(EIF)/n")

# ---------------------------------------------------------------------------
# 2. TMLE-minus-one-step diagnostic identities (same data, same learners)
# ---------------------------------------------------------------------------
cat("\n--- 2. TMLE vs one-step identities ---\n")

tm_lin <- modulatR::run_tmle_for_LMTP(
  ds1, make_factory(ds1),
  fluctuation = modulatR::LMTPFluctuationSubmodel$new(
    family = gaussian(), method = "linear"
  )
)
tm_wt <- modulatR::run_tmle_for_LMTP(
  ds1, make_factory(ds1),
  fluctuation = modulatR::LMTPFluctuationSubmodel$new(
    family = gaussian(), method = "weighted"
  )
)

# Nuisances identical across the three fits (deterministic GLM learners)
check(tm_lin$m_init$m_obs[[1]] %~% os$m_init$m_obs[[1]] &&
      tm_lin$omega[[1]] %~% os$omega[[1]],
      "TMLE and one-step share identical nuisances")

# 2a. Linear TMLE - one-step: the balance-factor identity.
# At tau = 1: TMLE = Pn[m_d] + eps * Pn[H_d] and onestep = Pn[m_d] + C,
# with C = Pn[omega (Y - m)] and eps = C / Pn[omega^2]. Since
# eps * Pn[H_d] = Pn[m_star_d] - Pn[m_init_d], the difference is computable
# entirely from stored fit slots:
omega <- os$omega[[1]]
m_obs <- os$m_init$m_obs[[1]]
C_hat <- mean(omega * (df1$Y - m_obs))

eps_lin <- unname(tm_lin$eps[[1]])
check(eps_lin %~% (C_hat / mean(omega^2)),
      "linear eps equals Pn[omega(Y - m)] / Pn[omega^2]")
check((tm_lin$estimate - os$estimate) %~%
        ((mean(tm_lin$m_star$m_d[[1]]) - mean(tm_lin$m_init$m_d[[1]])) - C_hat),
      "linear TMLE - one-step equals eps * Pn[H_d] - Pn[omega(Y - m)]")

# 2b. Weighted TMLE - one-step = C * (1/Pn[omega] - 1)
eps_wt <- unname(tm_wt$eps[[1]])
check(eps_wt %~% (C_hat / mean(omega)),
      "weighted eps equals Pn[omega(Y - m)] / Pn[omega]")
check((tm_wt$estimate - os$estimate) %~% (eps_wt - C_hat),
      "weighted TMLE - one-step equals eps - Pn[omega(Y - m)]")

# ---------------------------------------------------------------------------
# 3. Subgroup one-step
# ---------------------------------------------------------------------------
cat("\n--- 3. run_subgroup_onestep_for_LMTP ---\n")

ds <- make_ds(df)
subgroup_funs <- list(
  V0 = function(d) as.numeric(d$V == 0),
  V1 = function(d) as.numeric(d$V == 1)
)

os_sub <- modulatR::run_subgroup_onestep_for_LMTP(
  ds = ds,
  nuisance_factory = make_factory(ds),
  subgroup_funs = subgroup_funs
)

check(identical(names(os_sub$estimate), c("V0", "V1")),
      "subgroup estimates named")

# 3a. Estimate = group mean of phi
phi <- os_sub$omega[[1]] * (df$Y - os_sub$m_init$m_obs[[1]]) +
  os_sub$m_init$m_d[[1]]
check(unname(os_sub$estimate[["V1"]]) %~% mean(phi[df$V == 1]),
      "subgroup one-step equals within-group mean of phi")

# 3b. Per-subgroup scores solved
check(all(abs(colMeans(os_sub$eif)) < 1e-10),
      "Pn[D_j] = 0 for every subgroup")

# 3c. Overlapping subgroups accepted
os_ovl <- modulatR::run_subgroup_onestep_for_LMTP(
  ds = ds,
  nuisance_factory = make_factory(ds),
  subgroup_funs = list(
    V1     = function(d) as.numeric(d$V == 1),
    high_L = function(d) as.numeric(d$L > 0)
  )
)
check(all(abs(colMeans(os_ovl$eif)) < 1e-10),
      "overlapping subgroups: both scores solved")

# 3d. Consistency with the scalar runner on the V1 stratum:
# subgroup one-step for V1 uses pooled nuisances, so it will differ from the
# stratified scalar one-step; here we only check both are finite and close
# in this easy DGP.
check(is.finite(os_sub$estimate[["V1"]]) &&
        abs(os_sub$estimate[["V1"]] - os$estimate) < 1,
      "pooled-subgroup and stratified one-step agree to first order")

cat("\nAll tests passed.\n")
