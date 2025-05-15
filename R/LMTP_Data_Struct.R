
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
#'   data = df,
#'   id_col = 'id',
#'   n_timesteps = 2,
#'   exposure_cols = c('At1', 'At2'),
#'   time_varying_covariate_cols = list(c('L1t1', 'L2t1'), c('L1t2', 'L2t2')),
#'   outcome_col = 'Y'
#'   )
#'
LMTP_Data_Struct <- R6Class(
  "LMTP_Data_Struct",
  public = list(
    data = NULL,
    id_col = NULL,
    n_timesteps = NULL,
    exposure_cols = NULL,
    baseline_covariate_cols = NULL,
    time_varying_covariate_cols = NULL,
    outcome_col = NULL,
    Y = NULL,

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
      self$Y <- data[[outcome_col]]

      private$validate_input_data_structure()
    },

    H_cols = function(t) {
      if (! (is.numeric(t) && t %% 1 == 0) || length(t) != 1) {
        stop("t must be a length 1 integer.")
      }

      relevant_cols <-
        c(
          if (t == 1) NULL else self$exposure_cols[1:(t-1)],
          if (! is.null(self$baseline_covariate_cols)) self$baseline_covariate_cols[1:t] else NULL,
          if (! is.null(self$time_varying_covariate_cols)) {
            unlist(self$time_varying_covariate_cols[1:t])
          } else NULL
        )

      return(relevant_cols)
    },

    H = function(t) {
      self$data[,self$H_cols(t)]
    },

    A = function(t) {
      self$data[[self$exposure_cols[[t]]]]
    },

    L = function(t) {
      self$data[,self$time_varying_covariate_cols[[t]], drop=FALSE]
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

        if (! is.null(self$baseline_covariate_cols)) {
          if (! is.character(self$baseline_covariate_cols) || length(self$baseline_covariate_cols) != self$n_timesteps) {
            stop("baseline_covariate_cols must be a character vector of length n_timesteps")
          }
        }

        if (! is.null(self$time_varying_covariate_cols)) {
          if (! is.list(self$time_varying_covariate_cols) || length(self$time_varying_covariate_cols) != self$n_timesteps) {
            stop("time_varying_covariate_cols must be a list of length n_timesteps")
          }
        }

        if (! all(sapply(self$time_varying_covariate_cols, \(col_vec) {
          is.character(col_vec) }))) {
          stop("time_varying_covariate_cols must be a list of character vectors.")
        }

        if (! is.character(self$outcome_col) || length(self$outcome_col) != 1) {
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
