
# simulation details are passed in sbatch_99_run_everything.sh
commandLineArgs <- commandArgs()
task_i <- commandLineArgs[1]

library(here)

source(here("inst/simstudy1/code/00_setup.R"))
source(here("inst/simstudy1/code/01_simulate_LMTP.R"))
source(here("inst/simstudy1/code/03_sample_and_fit_lmtp.R"))
