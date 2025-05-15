#' lmtpiTask Class Constructor
#'
#' An lmtpiTask is an object that stores the relevant setup and methods for
#' estimating and constructing confidence intervals for a modified treatment
#' policy effect.
#'
#' @importFrom R6 R6Class
LMTPI_Task <- R6::R6Class("lmtpiTask",
  public = list(

    # initialize method for LMTPI_Task$new()
    initialize = function(
      data = NULL,
      policy = NULL,
      config = list(),
      network = NULL,
      exposure_map = NULL
    ) {
      self$data <- data
      self$policy <- policy
      self$config <- config
      self$network <- network
      self$exposure_map <- exposure_map
    },

    fit = function() {
      # ...
      # outcome_regression_model <- fit_outcome_regression_model()
      # treatment_model <- fit_treatment_model()
    }
  )
)
