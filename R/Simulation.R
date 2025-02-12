#' Simulation R6Class 
#' 
#' Takes a specification for a data generating process (`dgp`), 
#' list of `estimators`, `config`, list of `summary_statistics` to compute, and 
#' provides a `get_results` method. 
#'
#' The core idea is that a statistical simulation study consists of 
#' specifying a repeatable data generating process, some functions (estimators) 
#' to run on each generated data sample, and some summary statistics to compute 
#' from the simulation results (typically that indicate aspects of the performance
#' of the estimators considered). This is represented by the following pipeline:
#' 
#'   dgp          
#'   estimators    ->  Simulation$new( ... ) -> sim$run() -> sim$get_results() 
#'   config 
#'   summary_fns  
#' 
#' @examples 
#' # Example Usage
#' # Define a data generating process
#' dgp <- function(n) data.frame(x = rnorm(n), y = rnorm(n))
#' 
#' # Define some estimators
#' estimators <- list(
#'   mean_estimator = function(data) mean(data$x),
#'   var_estimator = function(data) var(data$x)
#' )
#' 
#' # Define a summary statistics function
#' summary_func <- function(iter = NULL, est_results, data = NULL) {
#'   data.frame(
#'     mean_est = est_results$mean_estimator,
#'     var_est = est_results$var_estimator
#'   )
#' }
#' 
#' # Create a simulation object
#' sim <- Simulation$new()
#' 
#' # Set up the simulation
#' sim$set_dgp(dgp)
#' sim$set_estimators(estimators)
#' sim$set_config(list(replications = 500, sample_size = 50))
#' sim$set_summary_stats(summary_func)
#' 
#' # Run the simulation
#' sim$run()
#' 
#' # Retrieve results
#' results <- sim$get_results()
#' head(results)
#' 
Simulation <- R6Class("Simulation",
  public = list(
    dgp = NULL,
    estimators = NULL,
    config = list(replications = 100, sample_size = 100),
    summary_stats = NULL,
    results = NULL,
    
    # Method to initialize the simulation object
    initialize = function() {},
    
    # Method to set the data generating process
    set_dgp = function(dgp_func) {
      if (!is.function(dgp_func)) {
        stop("dgp must be a function.")
      }
      self$dgp <- dgp_func
    },
    
    # Method to set the estimators
    set_estimators = function(estimator_list) {
      if (!is.list(estimator_list) || !all(sapply(estimator_list, is.function))) {
        stop("estimators must be a list of functions.")
      }
      self$estimators <- estimator_list
    },
    
    # Method to set the configuration
    set_config = function(config_list) {
      if (!is.list(config_list)) {
        stop("config must be a list.")
      }
      self$config <- modifyList(self$config, config_list)
    },
    
    # Method to set summary statistics
    set_summary_stats = function(summary_func) {
      if (!is.function(summary_func)) {
        stop("summary_stats must be a function.")
      }
      self$summary_stats <- summary_func
    },
    
    # Method to run the simulation
    run = function() {
      if (is.null(self$dgp) || is.null(self$estimators) || is.null(self$summary_stats)) {
        stop("Please set dgp, estimators, and summary_stats before running the simulation.")
      }
      
      if (! self$config$quiet) {
        message("Running simulation...")
      }

      replications <- self$config$replications
      sample_size <- self$config$sample_size
      sim_results <- vector("list", replications)
      
      for (i in seq_len(replications)) {
        data <- self$dgp(sample_size)
        est_results <- lapply(self$estimators, function(estimator) estimator(data))
        sim_results[[i]] <- self$summary_stats(i, est_results, data)
      }
      
      self$results <- do.call(rbind, sim_results)
    },
    
    # Method to retrieve results
    get_results = function() {
      if (is.null(self$results)) {
        stop("No results available. Run the simulation first.")
      }
      return(self$results)
    }
  )
)