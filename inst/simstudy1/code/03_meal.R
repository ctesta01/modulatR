library(here)

# the purpose of a meal is to digest its ingredients.
#
# we will consume the fit results from each of our jobs in a job array,
# and then run evaluation and visualization.

# dependencies
basedir <- here("inst/simstudy1/")
source(here(basedir, "code/01_run_on_cluster/00_dependencies.R"))

# assign job_id to be an NA_integer_
job_id <- NA_integer_
# define the sim_array_config
source(here(basedir, "code/01_run_on_cluster/01_job_array_config.R"))

# read in fit results
result_dirs <- list.files(here(basedir, "results/sim_results/"))
fit_results <- lapply(result_dirs, \(job_dir) {
  readRDS(here(basedir, 'results/sim_results', job_dir, 'fit_results.rds'))
})

# extract fit_results -- and expand into one giant dataframe
fit_results <- purrr::map_dfr(
  fit_results,
  \(x) {
    # append the job_id from each fit_results
    bind_cols(tibble(
      job_id     = x$job_id
    ),
    as_tibble(x$fit_results))
  }
)

fit_results$.rep <- as.integer(fit_results$.rep)

# bring in our setup data describing each job
# sim_array_config has important information describing every job.
#
# we left_join our fit_results into it by job_id so that we have
# that information available to us in evaluation and visualization.
fit_results <- dplyr::left_join(sim_array_config, fit_results, by = c('job_id' = 'job_id'))

# construct experiment ----------------------------------------------------

source(here(basedir, "code/01_run_on_cluster/02_dgp.R"))
source(here(basedir, "code/01_run_on_cluster/03_methods.R"))

# NOTE that the experiment defined here may differ from the one defined in
# R/01_run_on_cluster/04_experiment.R
#
# The reason for that is that the experiment that has run_experiment() called
# on it on the cluster may be just a subset of the whole experiment.
#
# The experiment we define here should include all of the dgps and methods
# used.
#
# This experiment is for:
#  - running the evaluators and visualizers, and
#  - rendering the documentation
#
experiment <- create_experiment(
  name = "LMTP Overall Effect") |>
  add_dgp(lmtp_dgp) |>
  add_method(m_modulatR)

# load evaluators and visualizers -----------------------------------------

source(here(basedir, "code/02_evaluators_and_visualizers/05_evaluators.R"))
source(here(basedir, "code/02_evaluators_and_visualizers/06_visualizers.R"))

# add evaluators and visualizers to recipe --------------------------------

# add the evaluator and visualizations to the experiment
experiment <- experiment |>
  add_evaluator(eval_lmtp) |>
  add_visualizer(viz_coverage_by_n) |>
  add_visualizer(viz_bias_by_n) |>
  add_visualizer(viz_sqrtn_error_hist)

# we have to call `save = TRUE` in order for documentation to render properly
experiment$evaluate(fit_results, save = TRUE)
experiment$visualize(fit_results, save = TRUE)

# setup and render documentation ------------------------------------------

# render the experiment's documentation
render_docs(experiment)

# View the rendered file now!
browseURL(here("results/LMTP methods/LMTP methods.html"))
