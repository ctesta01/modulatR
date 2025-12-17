
## Perform a forward simulation to get the true difference in counterfactual means.
##
## Only needs to be run once.
##
## Default is to use 2 million sample size for calculation.

compare_lmtp_to_natural <- function(n = 2e6L) {

  # Natural world (no policy)
  dat_natural <- simulate_lmtp_data_under_counterfactual_policy(
    n = n_truth, policy_seq = NULL)
  mean_Y_natural <- mean(dat_natural$Y)
  mean_Y_natural_V1 <- dat_natural |> dplyr::filter(V == 1) |> dplyr::pull(Y) |> mean()
  mean_Y_natural_V2 <- dat_natural |> dplyr::filter(V == 2) |> dplyr::pull(Y) |> mean()

  # LMTP world with shift policy
  shift_mtp   <- mtp_additive_shift(delta = -0.05)
  policy_seq  <- repeat_policy_over_time(shift_mtp, 3)

  dat_policy  <- simulate_lmtp_data_under_counterfactual_policy(
    n = n_truth, policy_seq = policy_seq)
  mean_Y_policy <- mean(dat_policy$Y)
  mean_Y_policy_V1 <- dat_policy |> dplyr::filter(V == 1) |> dplyr::pull(Y) |> mean()
  mean_Y_policy_V2 <- dat_policy |> dplyr::filter(V == 2) |> dplyr::pull(Y) |> mean()

  # Compare the results: Natural vs. Under Policy
  simulated_truth <- c(
    EY_natural = mean_Y_natural,
    EY_policy  = mean_Y_policy,
    EY_natural_V1 = mean_Y_natural_V1,
    EY_policy_V1 = mean_Y_policy_V1,
    EY_natural_V2 = mean_Y_natural_V2,
    EY_policy_V2 = mean_Y_policy_V2,
    diff       = mean_Y_policy - mean_Y_natural,
    diff_V1 = mean_Y_policy_V1 - mean_Y_natural_V1,
    diff_V2 = mean_Y_policy_V2 - mean_Y_natural_V2
  )
  simulated_truth

}

# run the calculation
true_counterfactual_means <- compare_lmtp_to_natural()

# save the results
writeRDS(true_counterfactual_means, here(prefix_dir, results_dir, "true_counterfactual_means.RDS"))
