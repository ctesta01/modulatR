# `modulatR` — Longitudinal Modified Treatment Policies Modulated by Subgroup Effects

`{modulatR}` aims to provide
  
  * a clean interface for specifying (longitudinal) modified treatment policies
  * with support for heterogeneous treatment effect estimation for (categorical/discrete) subgroups
  * and supporting functions that visualize / report on the estimates


## installation 

`{modulatR}` can be installed via GitHub using 

```r
devtools::install_github("ctesta01/modulatR")
```

## simple demos

```r
# redraft as of Apr 23 2025
# 
# after reconsideration, I think modularizing the interface further would be 
# an improvement. 

# workflow: 
# users should first specify the following: 

lmtp_data_struct <- ... 

config <- ... 

policy <- ... 

# and then users can run: 

lmtp_task$new(
  lmtp_data_struct,
  network,
  config,
  exposure_map,
  estimation_strategy = 'tmle'
  )

lmtp_task$run()
```


```r
library(modulatR)

# specify a modified treatment policy
mtp <- MTP$new(
    policy = function(a, covariates) { a - 5 }, 
    inverse_policy = function(A, L) { A - 5 },
    derivative_of_inverse_policy = function(A, L) { 1 }
    # default is to assume domain = c(-Inf, Inf)
    # the domain can be a collection (list) of subsets of the support of (A, L) on which 
    # the policy is smooth and invertible
)

load(lmtpi_dataset, package = 'lmtpi') # load example dataset
load(lmtpi_interference_matrix, package = 'lmtpi') # load example interference structure 

# To fit superlearned outcome models, the user specifies a list of 
# functions that accept a model specification, train that model, and then 
# return functions that predict from that model. 
# 
outcome_regression_models <- list(
  glm = function(X, Y) {
    model <- glm(Y ~ X)
    return(function(X_new) { predict(model, X_new) })},
  glmm = function(X, Y) {
    model <- lmer(Y ~ X + (id | group))
    return(function(X_new) { predict(model, X_new) } ) },
  rf = function(X, Y) {
    model <- randomForest(Y ~ X),
    return(function(X_new) { predict(model, X_new) }) }
)



# run the lmtpi algorithm 
task <- lmtpi::lmtpi_task(
    data = lmtpi_dataset,
    interference_matrix = lmtpi_interference_matrix,
    mtp = mtp,
    outcome_var = "cases", 
    exposure_var = "treatment",
    time_var = "year",
    confounders = c("confounder1", "confounder2", "confounder3"),
    estimation_config = 
      list(
        n_folds = 10,
        outcome_type = 'count',
        SL_library = list(
          # outcome_model = c('glm', 'glmm', 'randomForest', 'glmnet' ...),
          # propensity_ratio_model = c('glm', 'haldensify'))
          outcome_model = list(
            list(
              model = 'glm',
              formula = Y ~ ...,
              hyperparameters = ...
            ),
            list(
              model = 'lmer',
              formula = Y ~ (id | group) + ... ,
              hyperparameters = ...
            ),
            # closure of a model; 
            # just a function of X and Y 
          )
      )
    )
  
task$fit()
#> estimate: XX.XX
#> confidence intervals: [XX.XX, XX.XX] 
```

# Development Goals 

  - Write an MTP class satisfying the above specification 
  - Write vignettes that:
    - Show the most basic use 
  - Build out support for longitudinal estimation
  - Build out support for piecewise smooth invertible policies
  - Build out support for estimation with interference 
  - Build support for outcome types: continuous, binary, count 
  
# References and Similar Works

  - 
