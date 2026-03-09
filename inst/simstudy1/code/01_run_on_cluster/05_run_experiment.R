# R/05_run_experiment.R

library(here)

# setup specific to a job in a job-array
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1])

basedir <- here("inst/simstudy1/code/01_run_on_cluster")

source(here(basedir, "00_dependencies.R"))

# set up parallelization
#
# All of the job-specific configuration needs to be done _after_ setting
# up parallelization.
n_workers <- availableCores()

message(paste0("available cores detected: ", n_workers))

switch(Sys.info()[['sysname']],
       Windows = { plan(multisession, workers = n_workers) },
       Darwin = { plan(multisession, workers = n_workers) },
       Linux = { plan(multicore, workers = n_workers) }
)

# set random seed specific to each job
# seed_base <- 10000
# job_seed <- seed_base + job_id
# set.seed(job_seed)

source(here(basedir, "01_job_array_config.R"))
source(here(basedir, "02_dgp.R"))
source(here(basedir, "03_methods.R"))
source(here(basedir, "04_experiment.R"))

# Each task writes into its own directory to avoid collisions
out_dir <- file.path(here("inst/simstudy1/results"), "sim_results", sprintf("job_%04d", job_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# our experiment is defined in 04_experiment.R
fit_results <- fit_experiment(experiment, n_reps = n_reps_job,
                              future.globals = c('bind_rows'))

# save our simulation results
saveRDS(
  list(
    job_id = job_id,
    n_reps_job = n_reps_job,
    fit_results = fit_results
  ),
  file = file.path(out_dir, "fit_results.rds")
)

# this is where cluster work can end
