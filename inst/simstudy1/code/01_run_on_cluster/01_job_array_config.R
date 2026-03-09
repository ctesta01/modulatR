# R/01_sim_array_config.R

# In this file we define a config data.frame that controls for each
# simulation in the SLURM job-array, what should that job do.
#
# The details that I anticipate we may want to control include:
#
#   - What DGP to use
#   - What sample size to use
#   - How many replications to run
#
# While I am going to focus on these in this example, this script
# should be able to be extended to other config options in a similar manner.


# !!!!!!!
# you must ensure this matches the job-array given in slurm/slurm_run_experiment.sh
n_jobs <- 100
# I recommend that you calculate this number yourself by accounting for
#   - how many DGPs you have
#   - how many sample sizes you will be using
#   - for each combination of the above, how many replications will you run in each job

# here we will use 1 DGPs with 4 sample sizes and each job will run 100 replications.

sim_array_config <- data.frame(
  job_id = 1:n_jobs,
  dgp = rep(c('lmtp_dgp'), each = 100),
  sample_size = rep(c(300, 500, 1000, 10000), 25),
  n_reps_job = rep(10, 100))

# in the context of R/05_run_experiment.R, this script is being run for a
# specific simulation --- now get the details of that specific simulation setup.
dgp_choice <- sim_array_config[job_id, 'dgp']
sample_size <- sim_array_config[job_id, 'sample_size']
n_reps_job <- sim_array_config[job_id, 'n_reps_job']


# global simulation config parameters
delta <- -0.05         # policy shift in your example
sd_for_noise <- 10     # for specifying the noise level in simulations
