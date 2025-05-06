
#' Structured Data Object for Estimating Longitudinal Modified Treatment Policies
#'
#' The goal of having \code{LMTP_Data_Struct} is to have a standardized way to
#' refer to the time-varying covariates, exposure, and outcome data used in
#' estimating an LMTP causal parameter.
#'
#' @importFrom dplyr select mutate contains
#' @export
#'
#' @examples
#'
#' df <- data.frame(
#'   id = 1:5,
#'   At1 = rnorm(n=5),
#'   At2 = rnorm(n=5),
#'   L1t1 = rnorm(n=5),
#'   L2t1 = rnorm(n=5),
#'   L1t2 = rnorm(n=5),
#'   L2t2 = rnorm(n=5),
#'   Y = rnorm(n=5)
#'   )
#'
#' lmtp_data <- LMTP_Data_Struct$new(
#'   id_colname = 'id',
#'   n_timesteps = 2,
#'   exposures = df |> select(id, At1, At2),
#'   baseline_covariates = NULL,
#'   time_varying_covariates = lapply(1:2, \(t) { df |> select(id, contains(paste0('t', t))) }),
#'   outcomes = df |> select(id, Y)
#'   )
#'
#
#  LMTP_Data_Struct <- R6Class(
#    "LMTP_Data_Struct",
#    public = list(
#      id_colname = NULL,
#      n_timesteps = NULL,
#      exposures = NULL,
#      baseline_covariates = NULL,
#      time_varying_covariates = NULL,
#      outcomes = NULL,
#
#      initialize = function(
#        id_colname,
#        n_timesteps,
#        exposures,
#        time_varying_covariates,
#        baseline_covariates = NULL,
#        outcomes) {
#
#        self$id_colname <- id_colname
#        self$n_timesteps <- n_timesteps
#        self$exposures <- exposures
#        self$time_varying_covariates <- time_varying_covariates
#        self$baseline_covariates <- baseline_covariates
#        self$outcomes <- outcomes
#
#        private$validate()
#      }
#    ),
#
#    private = list(
#      validate_that_id_colname_appears_everywhere = function() {
#        # check the id_colname appears in the exposures data given
#        if (! is.list(self$exposures) ||
#          ! all(sapply(self$exposures, \(df) self$id_colname %in% colnames(df))) {
#          stop("the provided id_colname does not appear in the exposures colnames.")
#        }
#
#        # check the id_colname appears in the baseline covariates data given
#        if (! is.null(self$baseline_covariates)) {
#          if (! self$id_colname %in% colnames(self$baseline_covariates)) {
#            stop("the provided id_colname does not appear in the baseline_covariates colnames.")
#          }
#        }
#
#        # check that the time_varying_covariates are given as a list of data frames,
#        # one for each timestep.
#        if (! is.list(self$time_varying_covariates) ||
#            length(self$time_varying_covariates) != self$n_timesteps ||
#            ! all(sapply(self$time_varying_covariates, is.data.frame))) {
#          stop("The time_varying_covariates must be entered as a list of dataframes,
#  one for each time-step.")
#        }
#
#        # check that the time varying covariates have the id_colname
#        if (! all(sapply(self$time_varying_covariates,
#                         \(df) { self$id_colname %in% colnames(df) }))
#        ) {
#          stop("The id_colname does not appear in every one of the time_varying_covariates dataframes provided.")
#        }
#
#        # check that id_colname appears in the outcomes
#        if (! self$id_colname %in% colnames(self$outcomes)) {
#          stop("The id_colname does not appear in every one of the time_varying_covariates dataframes provided.")
#        }
#
#        return(invisible(TRUE))
#      },
#      validate_data_are_ordered_properly = function() {
#
#        # check the data are ordered properly
#        id_ordering <- self$exposures[[self$id_colname]]
#        num_unique_ids <- length(unique(id_ordering))
#        if (num_unique_ids != length(id_ordering)) {
#          stop("The ids in the id_column given must be unique. This is not the case for the exposure data.")
#        }
#
#        if (! is.null(self$baseline_covariates)) {
#          if (num_unique_ids != length(self$baseline_covariates[[self$id_colname]])) {
#            stop("The ids in the id_column given must be unique and match across the
#    exposures, baseline_covariates, time_varying_covariates, and outcomes.
#    This is not the case for the baseline_covariates.")
#          }
#        }
#
#        if (any(sapply(self$time_varying_covariates, \(df) num_unique_ids != length(df[[self$id_colname]])))) {
#          stop("The ids in the id_column given must be unique and match across the
#  exposures, baseline_covariates, time_varying_covariates, and outcomes.
#  This is not the case for the time_varying_covariates.")
#        }
#
#        if (num_unique_ids != length(self$outcomes[[self$id_colname]])) {
#          stop("The ids in the id_column given must be unique and match across the
#  exposures, baseline_covariates, time_varying_covariates, and outcomes.
#  This is not the case for the outcomes.")
#        }
#
#        # make sure baseline covariates are ordered properly
#        if (! is.null(self$baseline_covariates)) {
#          if (! all(id_ordering == self$baseline_covariates[[self$id_colname]])) {
#            self$baseline_covariates[[id_colname]] <-
#              factor(self$baseline_covariates, levels = id_ordering)
#            self$baseline_covariates %<>% arrange({{ id_colname}})
#          }
#        }
#
#        # make sure outcomes are ordered properly
#        if (! all(id_ordering == self$outcomes[[self$id_colname]])) {
#          self$outcomes[[id_colname]] <-
#            factor(self$outcomes, levels = id_ordering)
#          self$outcomes %<>% arrange({{ id_colname}})
#        }
#
#        # make sure time-varying-covariates are ordered properly
#        for (t in 1:self$n_timesteps) {
#          if (! all(id_ordering == self$time_varying_covariates[[t]][[self$id_colname]])) {
#            self$time_varying_covariates[[t]][[id_colname]] <-
#              factor(self$time_varying_covariates[[t]], levels = id_ordering)
#            self$time_varying_covariates[[t]] %<>% arrange({{ id_colname}})
#          }
#        }
#
#        return(invisible(NULL))
#      },
#      validate = function() {
#        private$validate_that_id_colname_appears_everywhere()
#        private$validate_data_are_ordered_properly()
#        return(invisible(NULL))
#      }
#    )
#  )
LMTP_Data_Struct <- R6Class(
  "LMTP_Data_Struct",
  public = list(
    data = NULL,
    id_col = NULL,
    n_timesteps = NULL,
    exposure_cols = NULL,
    baseline_covariate_cols = NULL,
    time_varying_covariate_cols = NULL,
    outcome_cols = NULL,

    initialize = function(
      data,
      id_col,
      n_timesteps,
      exposure_cols,
      baseline_covariate_cols = NULL,
      time_varying_covariate_cols = NULL,
      outcome_col) {

      # throughout, col refers to the name of one column;
      # cols refers to the names of multiple columns
      self$data <- data
      self$id_col <- id_col
      self$n_timesteps <- n_timesteps
      self$exposure_cols <- exposure_cols
      self$baseline_covariate_cols <- baseline_covariate_cols
      self$time_varying_covariate_cols <- time_varying_covariate_cols
      self$outcome_col <- outcome_col

      validate_input_data_structure()
    }
  ),
  private = list(
    validate_input_data_structure = function() {
      if (! is.numeric(self$n_timesteps) && self$n_timesteps %% 1 == 0) {
        stop("n_timesteps must be an integer.")
      }

      if (! self$id_col %in% colnames(self$data)) {
        stop("id_col must appear among the column names of the data given.")
      }

      # TODO:
      # write code that standardizes structure for shortened input format for
      # single-time-point settings


      # check that the dimensions of exposure_cols, baseline_covariate_cols,
      # time_varying_covariate_cols, and outcome_col are correct
      # =====

      if (self$n_timesteps > 1) {
        if (! is.character(self$exposure_cols) || length(self$exposure_cols) != self$n_timesteps) {
          stop("exposure_cols must be a character vector of length n_timesteps")
        }

        if (! is.character(self$baseline_covariate_cols) || length(self$baseline_covariate_cols) != self$n_timesteps) {
          stop("baseline_covariate_cols must be a character vector of length n_timesteps")
        }

        if (! is.list(self$time_varying_covariate_cols) || length(self$baseline_covariate_cols) != self$n_timesteps) {
          stop("time_varying_covariate_cols must be a list of length n_timesteps")
        }

        if (! all(sapply(self$time_varying_covariate_cols, \(col_vec) {
          is.character(col_vec) }))) {
          stop("time_varying_covariate_cols must be a list of character vectors.")
        }

        if (! is.character(self$outcome_col) || length(self$outcome_col) == 1) {
          stop("outcome_col must be a character vector of length 1.")
        }


        # check that the exposure_cols, baseline_covariate_cols, time_varying_covariate_cols,
        # outcome_col appear in the data
        # =====

        # helper function to check that a character vector appears as colnames
        # of the data
        appears_in_data_chr_vec_input <- function(chr_vec) {
          return(all(chr_vec %in% colnames(self$data)))
        }

        # helper function that a list of character vector appears as colnames of the data
        appears_in_data_list_of_chr_vec_input <- function(list_of_chr_vec) {
          return(all(sapply(list_of_chr_vec, \(chr_vec) {
            appears_in_data_chr_vec_input(chr_vec)
          })))
        }

        if (! appears_in_data_chr_vec_input(self$outcome_col)) {
          stop("the outcome_col must appear in the data colnames")
        }
        if (! appears_in_data_list_of_chr_vec_input(self$exposure_cols)) {
          stop("the outcome_cols must appear in the data colnames")
        }
        if (! appears_in_data_list_of_chr_vec_input(self$baseline_covariate_cols)) {
          stop("the baseline_covariate_cols must appear in the data colnames")
        }
        if (! appears_in_data_list_of_chr_vec_input(self$time_varying_covariate_cols)) {
          stop("the time_varying_covariate_cols must appear in the data colnames")
        }
      }

      return(invisible(NULL))
    }
  )
)
