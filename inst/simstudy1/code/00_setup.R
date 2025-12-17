library(here)
library(stringr)
library(lubridate)
library(dplyr)
devtools::load_all(".")

prefix_dir <- here("inst/simstudy1")
date_str <- lubridate::today()
results_str <- str_c("results-", date_str)
dir.create(here(prefix_dir, results_str), showWarnings = FALSE)
dir.create(here(prefix_dir, str_c(results_str, "/sims")), showWarnings = FALSE)
dir.create(here(prefix_dir, str_c(results_str, "/agg")), showWarnings = FALSE)
dir.create(here(prefix_dir, str_c(results_str, "/sbatch-out")), showWarnings = FALSE)

# --- configuration ---
n_vals <- c(100, 300, 500, 1000, 10000)
n_reps_total <- 1000   # total repetitions per n
delta <- -0.05         # policy shift in your example
sd_for_noise <- 10     # for specifying the noise level in simulations

# simulations to be run
n_vals_per_sim <- rep(n_vals, each = n_reps_total)

# chunk the simulations into 50 tasks
# we know this is an even division, so no need for float math
n_sims_per_task <- length(n_vals_per_sim) / 50

sims_per_task <- lapply(1:50, function(task_i) {
  n_sims_per_task * (task_i - 1) + seq(1:n_sims_per_task)
})

# simulation settings for using randomForest
options(future.globals.maxSize = Inf) # if using lnr_rf_binary
