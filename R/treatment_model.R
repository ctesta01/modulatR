
#' @param density_ratio_approach One of either 'naive' or 'classification'
modulatR_ipw_estimator <- function(
    lmtp_data,
    policy,
    learners = list(nadir::lnr_glm),
    extra_super_learner_arguments = NULL,
    density_ratio_approach = 'classification'
) {

  n_timesteps <- lmtp_data$n_timesteps

  if (density_ratio_approach == 'naive') {

    density_ratio_denominators <- rep(1, nrow(lmtp_data$data))
    density_ratio_numerators <- rep(1, nrow(lmtp_data$data))

    for (timestep_t in n_timesteps:1) {

      A_col <- lmtp_data$exposure_cols[[timestep_t]]

      training_dataset <-
        dplyr::bind_cols(
          A = lmtp_data$A(timestep_t),
          Ht)

      # regress Y on everything else in the training data
      model_t <-
        do.call(
          what = nadir::super_learner,
          args = c(list(
            data = training_dataset,
            formulas = as.formula(paste0(A_col, " ~ .")),
            learners = learners,
            y_variable = A_col),
            extra_super_learner_arguments
          )
        )

      density_ratio_denominators <- density_ratio_denominators *
        model_t(lmtp_data$A(timestep_t))


    }
  } else if (density_ratio_approach == 'classification') {
  }

}
