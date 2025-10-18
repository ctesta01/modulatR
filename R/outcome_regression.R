
#'
modulatR_outcome_regression <- function(
    lmtp_data,
    policy,
    learners = list(nadir::lnr_glm),
    extra_super_learner_arguments = NULL
    ) {

  # get the name of the outcome variable
  Y_col <- lmtp_data$outcome_col

  # number of timesteps that covariates/exposure are observed
  n_timesteps <- lmtp_data$n_timesteps

  for (timestep_t in (n_timesteps):1) {

    # get history up to time t
    # Ht includes (L1, ..., Lt, A1, ..., A_{t-1})
    # We can think of Ht as the history up to (before, not including) At
    Ht = lmtp_data$H(timestep_t)
    Lt = lmtp_data$L(timestep_t)

    if (timestep_t == n_timesteps) { # for the first g-computation step

      # get the full history of exposures and covariates (At, Ht)
      training_dataset <-
        dplyr::bind_cols(
          Y = lmtp_data$Y,
          A = lmtp_data$A(timestep_t),
          Ht)


    } else if (timestep_t < n_timesteps & timestep_t >= 1) { # for subsequent g-computation steps

      # form a training dataset using the pseudo-outcome, At, Ht
      training_dataset <-
        dplyr::bind_cols(Y = pseudo_outcome,
                         A = lmtp_data$A(timestep_t),
                         Ht)
    }

    # make sure everything is named properly
    A_col = lmtp_data$exposure_cols[timestep_t]
    colnames(training_dataset)[1] <- Y_col
    colnames(training_dataset)[2] <- A_col

    # regress Y on everything else in the training data
    model_t <-
      do.call(
        what = nadir::super_learner,
        args = c(list(
          data = training_dataset,
          formulas = as.formula(paste0(Y_col, " ~ .")),
          learners = learners,
          y_variable = Y_col),
          extra_super_learner_arguments
          )
        )

    # use the MTP to apply the policy to A
    mtp_shifted_A_t <- policy$apply_policy(A = lmtp_data$A(timestep_t), L = lmtp_data$L(timestep_t))

    # combine into a dataframe of predictor variables to predict on
    prediction_dataset_under_policy_shift <- dplyr::bind_cols(A_col = mtp_shifted_A_t, Ht)

    # make sure the At variable is named properly
    colnames(prediction_dataset_under_policy_shift)[1] <- A_col

    # predict under the intervention at time t
    pseudo_outcome <- model_t(prediction_dataset_under_policy_shift)
  }

  return(pseudo_outcome)
}
