
# these are the sims task_i (i in the array of 1:50) is responsible for;
# this will be a vector of sample sizes to use.
sims_for_task_i <- sims_per_task[[task_i]]

set.seed(task_i) # this way every task in the job-array has a different seed

results <- list() # data struct to store results

for (j in 1:length(sims_for_task_i)) {
  # get sample size for jth simulation in task_i
  sample_n <- sims_for_task_i[[j]]

  # create data structure
  ds <- simulate_lmtp_data_under_counterfactual_policy(
    n = sample_n,
    policy_seq = NULL,
    sd_for_noise = sd_for_noise)

  # structure it in the LMTP Data Structure
  ds <- LMTP_Data_Struct$new(
    data = ds,
    id_col     = "id",
    n_timesteps = 3,
    A_cols     = c("A1", "A2", "A3"),
    L_cols     = list(c("L11", "L12"),
                      c("L21", "L22"),
                      c("L31", "L32")),
    W_cols     = c("W1", "V"),
    Y_col      = "Y")

  # define a sequence of policies that reduce exposure by -.05
  shift_mtp <- mtp_additive_shift(delta = delta)
  policy_seq <- repeat_policy_over_time(shift_mtp, 3)

  # specify a nuisance factory
  nuis <- LMTPNuisanceFactory$new(
    learners_g = nadir::lnr_lm_density,
    policy_seq = policy_seq,
    A_type = "continuous",
    repeat_fmls_lnrs_args = TRUE,
    g_mode = 'density'
  )

  # fit the model under the intervention policy for E[Y^d]
  fit_intv <- fit_tmle_for_LMTP(
    ds,
    policy_seq = policy_seq,
    learners_Q = nadir::lnr_lm,
    learners_g_factory = nuis,
    outcome_link = 'identity',
    repeat_lnrs = TRUE,
    method = 'tmle')

  # define the identity MTP and a sequence of identity MTPs
  identity_mtp <- mtp_additive_shift(delta = 0)
  policy_seq_identity <- repeat_policy_over_time(identity_mtp, 3)

  # fit the model under no intervention for E[Y]
  fit_natural <- fit_tmle_for_LMTP(
    ds,
    policy_seq = policy_seq_identity,
    learners_Q = nadir::lnr_lm,
    learners_g_factory = nuis,
    outcome_link = 'identity',
    repeat_lnrs = TRUE)

  # save the results
  results[[j]] <- list(
    intv_psi = fit_intv$psi,
    nat_psi = fit_natural$psi,
    intv_ci_low = fit_intv$ci95[1],
    intv_ci_high = fit_intv$ci95[2],
    nat_ci_low = fit_natural$ci95[1],
    nat_ci_high = fit_natural$ci95[2],
    diff_psi = fit_intv$psi - fit_natural$psi,
    diff_ci_low = quantile(fit_intv$ic - fit_natural$ic, 0.025, na.rm=TRUE),
    diff_ci_high = quantile(fit_intv$ic - fit_natural$ic, 0.975, na.rm=TRUE)
  )
}

# write the results to file
saveRDS(
  results,
  here(prefix_dir, results_str, "/sims/", "task_", task_i, "-sim_results.rds"))
