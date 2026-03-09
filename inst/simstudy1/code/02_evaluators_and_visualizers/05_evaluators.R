
if (! file.exists(
  here(basedir, 'results/true_counterfactual_means.RDS'))) {

compare_lmtp_to_natural <- function(n = 2e6L) {

  # Natural world (no policy)
  dat_natural <- lmtp_dgp$dgp_fun(
    n = n, policy_seq = NULL)$df
  mean_Y_natural <- mean(dat_natural$Y)
  mean_Y_natural_V1 <- dat_natural |> dplyr::filter(V == 1) |> dplyr::pull(Y) |> mean()
  mean_Y_natural_V2 <- dat_natural |> dplyr::filter(V == 2) |> dplyr::pull(Y) |> mean()

  # LMTP world with shift policy
  shift_mtp   <- mtp_additive_shift(delta = -0.05)
  policy_seq  <- repeat_policy_over_time(shift_mtp, 3)

  dat_policy  <- lmtp_dgp$dgp_fun(
    n = n, policy_seq = policy_seq)$df
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
saveRDS(true_counterfactual_means, here(basedir, "results/true_counterfactual_means.RDS"))

} else {
  true_counterfactual_means <- readRDS(here(basedir, "results/true_counterfactual_means.RDS"))
}





# define evaluators -------------------------------------------------------

true_effect <- true_counterfactual_means['EY_policy'] - true_counterfactual_means['EY_natural']



eval_lmtp <- simChef::create_evaluator(
  .eval_fun = function(fit_results, vary_params = NULL, ...) {
    fit_results |>
      group_by(sample_size) |>
      summarize(
        n_reps = dplyr::n(),
        coverage = mean(true_effect >= diff_ci_low & true_effect <= diff_ci_high, na.rm=TRUE),
        bias = mean(diff_psi - true_effect, na.rm=TRUE),
        sqrt_n_bias = mean(sqrt(sample_size) * (diff_psi - true_effect), na.rm=TRUE))
}, .name = 'lmtp_diff_eval')

