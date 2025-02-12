#' Super Learner: Cross-Validation Based Ensemble Learning
#'
#' @examples
#' learners <- list(
#'      glm = function(data, regression_formula = regression_formula) { model <- lm(formula = regression_formula, data = data); return(function(Xnew) { predict(model, newdata = Xnew) }) },
#'      rf = function(data, regression_formula) { model <- randomForest::randomForest(formula = regression_formula, data = data); return(function(Xnew) { predict(model, newdata = Xnew) }) },
#'      glmnet = function(data, regression_formula = regression_formula) {
#'        xvars <- attr(terms(regression_formula), "term.labels")
#'        yvar <- as.character(regression_formula)[2]
#'        model <- glmnet::glmnet(y = data[[yvar]], x = as.matrix(data[,xvars]), lambda = .2)
#'        return(function(Xnew) { as.vector(predict(model, newx = as.matrix(Xnew[,xvars]), type = 'response')) })
#'        }
#'   )
#'
#' # mtcars example ---
#' sl_model <- lmtp::super_learner(
#'   data = mtcars,
#'   regression_formula = mpg ~ cyl + hp,
#'   learners = learners)
#'
#' sl_model_predictions <- sl_model(mtcars)
#' fit_individual_learners <- lapply(learners, function(learner) { learner(data = mtcars, regression_formula = mpg ~ cyl + hp) } )
#' individual_learners_rmse <- lapply(fit_individual_learners, function(fit_learner) { rmse(fit_learner(mtcars) - mtcars$mpg) })
#'
#' print(paste0("super-learner rmse: ", rmse(sl_model_predictions - mtcars$mpg)))
#' individual_learners_rmse
#'
#' # iris example ---
#' sl_model <- super_learner(
#'   data = iris,
#'   regression_formula = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
#'   learners = learners)
#'
#' sl_model_predictions <- sl_model(iris)
#' fit_individual_learners <- lapply(learners, function(learner) { learner(data = iris, regression_formula = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width) } )
#' individual_learners_rmse <- lapply(fit_individual_learners, function(fit_learner) { rmse(fit_learner(iris) - iris$Sepal.Length) })
#'
#' print(paste0("super-learner rmse: ", rmse(sl_model_predictions - iris$Sepal.Length)))
#' individual_learners_rmse
#'
#' @export
super_learner <- function(
  data,
  learners,
  regression_formula,
  n_folds = 5,
  determine_superlearner_weights = determine_superlearner_weights_nnls,
  continuous_or_discrete = 'continuous') {

  # add a folds column
  data %<>% make_folds(n_folds = n_folds)

  # split into an n_folds length list of training_data and validation_data:
  # training_data contains n_folds-1 folds of data, validation_data contains 1 fold
  training_data <- lapply(
    1:n_folds, function(i) {data |> dplyr::filter(fold != i)})
  validation_data <- lapply(
    1:n_folds, function(i) {data |> dplyr::filter(fold == i)})

  # make a tibble/dataframe to hold the trained learners:
  # one for each combination of a specific fold and a specific model
  trained_learners <- tibble::tibble(
    split = rep(1:n_folds, length(learners)),
    learner_name = rep(names(learners), each = n_folds))

  # for each i in 1:n_folds and each model, train the model
  trained_learners$learned_predictor <- lapply(
    1:nrow(trained_learners), function(i) {
      learners[[trained_learners[[i,'learner_name']]]](
        data = training_data[[trained_learners[[i,'split']]]],
        regression_formula = regression_formula
      )
    }
  )

  # predict from each fold+model combination on the held-out data
  trained_learners$predictions_for_testset <- lapply(
    1:nrow(trained_learners), function(i) {
      trained_learners[[i,'learned_predictor']][[1]](validation_data[[trained_learners[[i, 'split']]]])
    }
  )

  # from here forward, we just need to use the split + model name + predictions on the test-set
  # to regress against the held-out (validation) data to determine the ensemble weights
  second_stage_SL_dataset <- trained_learners[,c('split', 'learner_name', 'predictions_for_testset')]

  # pivot it into a wider format, with one column per model, with columnname model_name
  second_stage_SL_dataset <- tidyr::pivot_wider(
    second_stage_SL_dataset,
    names_from = 'learner_name',
    values_from = 'predictions_for_testset')

  # Extract the Y-variable (its character name)
  #
  # This only supports simple Y variables, nothing like a survival right-hand-side or
  # a transformed right-hand-side.
  #
  # TODO: Add an error if Y is not a simple outcome
  formula_as_character <- as.character(regression_formula) # typically c("~", "yvar", "x_variables + ...")
  y_variable <- formula_as_character[[2]]

  # insert the validation data in another column next to the predictions
  # TODO: Change validation_data to Y
  #
  second_stage_SL_dataset$Y <- lapply(1:nrow(second_stage_SL_dataset), function(i) {
    validation_data[[second_stage_SL_dataset[[i, 'split']]]][[y_variable]]
  })

  # unnest all of the data (each cell prior to this contained a vector of either
  # predictions or the validation data)
  second_stage_SL_dataset <- tidyr::unnest(second_stage_SL_dataset, cols = everything())

  # drop the split column so we can simplify the following regression formula
  second_stage_SL_dataset$split <- NULL

  # regress the validation data on the predictions from every model with no intercept.
  # notice this is now for all of the folds
  #
  # TODO: Here we assume a continuous Y-variable and use a linear regression to
  # determine the SuperLearner weights;  we may want to support other types of
  # Y-variables like binary, count, and survival.  My theory on how to support
  # these most flexibly is to abstract the logic of the
  # model-weight-determination to a secondary function that eats
  # second_stage_SL_dataset and produces weights; that way the user can swap out
  # whatever they'd like instead, but several handy defaults are supported and
  # already coded up for users.
  #
  # TODO: An option for handling count outcomes / weighting the
  # outcomes/observations -- What may be a solution is multiplying the
  # rows by the square root of the desired weights...
  learner_weights <- determine_superlearner_weights(second_stage_SL_dataset)
  if (continuous_or_discrete == 'continuous') {
    # nothing needs to be done; leave the learner_weights as-is
  } else if (continuous_or_discrete == 'discrete') {
    max_learner_weight <- which(learner_weights == max(learner_weights))
    if (length(max_learner_weight) > 1) {
      warning("Multiple learners were tied for the maximum weight. Since discrete super-learner was specified, the first learner with the maximum weight will be used.")
      learner_weights <- rep(0, length(learner_weights))
      learner_weights[max_learner_weight[1]] <- 1
    }
  } else {
    stop("Argument continuous_or_discrete must be one of 'continuous' or 'discrete'")
  }

  # fit all of the learners on the entire dataset
  fit_learners <- lapply(
    1:length(learners), function(i) {
      learners[[i]](data = data, regression_formula = regression_formula)
    })

  # construct a function that predicts using all of the learners combined using
  # SuperLearned weights
  #
  # this is a closure that will be returned from this function
  predict_from_super_learned_model <- function(newdata) {
    # for each model, predict on the newdata and apply the model weights
    lapply(1:length(fit_learners), function(i) {
      fit_learners[[i]](newdata) * learner_weights[[i]]
    }) %>%
      Reduce(`+`, .) # aggregate across the weighted model predictions
  }

  return(predict_from_super_learned_model)
}


#' Determine SuperLearner Weights with Nonnegative Least Squares
#'
#' This function accepts a dataframe that is structured to have
#' one column `Y` and other columns with unique names corresponding to
#' different model predictions for `Y`, and it will use nonnegative
#' least squares to determine the weights to use for a SuperLearner.
#'
#' @export
determine_superlearner_weights_nnls <- function(data) {
  # use nonlinear least squares to produce a weighting scheme
  nnls_output <- nnls::nnls(
    A = as.matrix(dplyr::select(data, -Y)),
    b = data$Y)

  return(nnls_output$x)
}

#' Add a new column called folds to the data passed
#'
#' A helper function for inserting a column called `fold` into the
#' data with integers randomly sampled with replacement from `1:n_folds`.
#'
make_folds <- function(data, n_folds = 5) {
  if ('fold' %in% colnames(data)) {
    stop("The data passed to make_folds already has a folds column")
  }
  data$fold <- sample.int(n = n_folds, size = nrow(data), replace = TRUE)
  return(data)
}
