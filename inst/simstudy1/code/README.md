

The goal for this simulation study is to conduct a simple proof of concept that
the package is working correctly for one of its most basic applications:  estimating 
counterfactual means in the setting of observed time-varying confounding. 

```
code/
  00_setup.R                          # defines simulation structure -- needs to be run by every file 
  01_simulate_LMTP.R                  # used for all simulations
  02_get_true_counterfactual_means.R  # only used to determine true E[Y^d], E[Y] once 
  ...
  98_calculate_true_counterfactual_means_job.R  # a separate task to be run 
  99_run_everything.R    # a parameterized job (array_id = 1-50) to run 80 simulations of various sizes 
results/

```

