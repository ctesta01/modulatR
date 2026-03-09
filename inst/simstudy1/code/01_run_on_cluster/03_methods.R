
m_modulatR_fun <- function(df) {

    sample_n <- nrow(df)

    # structure it in the LMTP Data Structure
    ds <- modulatR::LMTP_Data_Struct$new(
      data = df,
      id_col     = "id",
      n_timesteps = 3,
      A_cols     = c("A1", "A2", "A3"),
      L_cols     = list(c("L11", "L12"),
                        c("L21", "L22"),
                        c("L31", "L32")),
      W_cols     = c("W1", "V"),
      Y_col      = "Y")

    # define a sequence of policies that reduce exposure by -.05
    shift_mtp <- modulatR::mtp_additive_shift(delta = -.05)
    policy_seq <- modulatR::repeat_policy_over_time(shift_mtp, 3)

    # specify a nuisance factory
    nuis <- LMTPNuisanceFactory$new(
      learners_g = nadir::lnr_lm_density,
      policy_seq = policy_seq,
      A_type = "continuous",
      repeat_fmls_lnrs_args = TRUE,
      g_mode = 'density'
    )

    # fit the model under the intervention policy for E[Y^d]
    fit_intv <- modulatR::run_tmle_for_LMTP(
      ds,
      policy_seq = policy_seq,
      learners_Q = nadir::lnr_lm,
      learners_g_factory = nuis,
      outcome_link = 'identity',
      repeat_lnrs = TRUE,
      method = 'tmle')

    # define the identity MTP and a sequence of identity MTPs
    identity_mtp <- modulatR::mtp_additive_shift(delta = 0)
    policy_seq_identity <- modulatR::repeat_policy_over_time(identity_mtp, 3)

    nuis_nat <- modulatR::LMTPNuisanceFactory$new(
      learners_g = nadir::lnr_lm_density,
      policy_seq = policy_seq_identity,
      A_type = "continuous",
      repeat_fmls_lnrs_args = TRUE
    )

    # fit the model under no intervention for E[Y]
    fit_natural <- modulatR::run_tmle_for_LMTP(
      ds,
      policy_seq = policy_seq_identity,
      learners_Q = nadir::lnr_lm,
      learners_g_factory = nuis_nat,
      outcome_link = 'identity',
      repeat_lnrs = TRUE)

    # save the results
    results <- list(
      intv_psi = fit_intv$psi,
      nat_psi = fit_natural$psi,
      intv_ci_low = fit_intv$ci95[1],
      intv_ci_high = fit_intv$ci95[2],
      nat_ci_low = fit_natural$ci95[1],
      nat_ci_high = fit_natural$ci95[2],
      diff_psi = fit_intv$psi - fit_natural$psi,
      diff_ci_low = (fit_intv$psi - fit_natural$psi) + qnorm(0.025)*sqrt(var(fit_intv$ic - fit_natural$ic)/sample_n),
      diff_ci_high = (fit_intv$psi - fit_natural$psi) + qnorm(0.975)*sqrt(var(fit_intv$ic - fit_natural$ic)/sample_n)
    )
  return(results)
}




# create simChef methods --------------------------------------------------

m_modulatR <- create_method(m_modulatR_fun,
                   .name = 'modulatR overall lmtp estimator')
