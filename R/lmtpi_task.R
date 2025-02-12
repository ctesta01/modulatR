#' lmtpiTask Class Constructor
#' 
#' An lmtpiTask is an object that stores the relevant setup and methods for 
#' estimating and constructing confidence intervals for a modified treatment 
#' policy effect. 
#' 
lmtpiTask <- R6Class("lmtpiTask",
  public = list(
    data = NULL,
    mtp = NULL,
    outcome_var = NULL,
    exposure_var = NULL,
    time_var = NULL,
    confounders = NULL,
    estimation_config = list(),

    # initialize method for lmtpTask$new() 
    initialize = function(
      data = NULL,
      mtp = NULL,
      outcome_var = NULL,
      exposure_var = NULL,
      time_var = NULL,
      confounders = NULL,
      estimation_config = list()
    ) {
      self$data <- data
      self$mtp <- mtp
      self$outcome_var <- outcome_var
      self$exposure_var <- exposure_var
      self$time_var <- time_var
      self$confounders <- confounders
      self$estimation_config <- estimation_config
    },

    fit = function() {
      # ... 
      outcome_regression_model <- fit_outcome_regression_model()
      treatment_model <- fit_treatment_model() 
    }
  )
)