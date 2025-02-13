#' Super Learner: Cross-Validation Based Ensemble Learning
#'
#' @examples
#' \dontrun{
#' library(lme4)
#' library(randomForest)
#'
#' learners <- list(
#'      glm = function(data, regression_formula, ...) {
#'        model <- lm(formula = regression_formula, data = data)
#'        return(function(newdata) { predict(model, newdata = newdata) })
#'        },
#'      rf = function(data, regression_formula, ...) {
#'        model <- randomForest::randomForest(formula = regression_formula, data = data)
#'        return(function(newdata) { predict(model, newdata = newdata) })
#'        },
#'      glmnet = function(data, regression_formula, ...) {
#'        xvars <- attr(terms(regression_formula), "term.labels")
#'        yvar <- as.character(regression_formula)[2]
#'        model <- glmnet::glmnet(y = data[[yvar]], x = as.matrix(data[,xvars]), lambda = .2)
#'        return(function(newdata) { as.vector(predict(model, newx = as.matrix(newdata[,xvars]), type = 'response')) })
#'        },
#'      lme4 = function(data, regression_formula, ...) {
#'        model <- lme4::lmer(formula = regression_formula, data = data)
#'        return(function(newdata) { predict(model, newdata = newdata) })
#'      }
#'   )
#'
#' # mtcars example ---
#' regression_formulas <- c(
#'   rep(c(mpg ~ cyl + hp), 3), # first three models use same formula
#'   mpg ~ (1 | cyl) + hp # lme4 uses different language features
#'   )
#'
#' sl_model <- super_learner(
#'   data = mtcars,
#'   regression_formula = regression_formulas,
#'   learners = learners)
#'
#' sl_model_predictions <- sl_model(mtcars)
#' fit_individual_learners <- lapply(1:length(learners), function(i) { learners[[i]](data = mtcars, regression_formula = regression_formulas[[i]]) } )
#' individual_learners_rmse <- lapply(fit_individual_learners, function(fit_learner) { rmse(fit_learner(mtcars) - mtcars$mpg) })
#'
#' print(paste0("super-learner rmse: ", rmse(sl_model_predictions - mtcars$mpg)))
#' individual_learners_rmse
#'
#' # iris example ---
#' sl_model <- super_learner(
#'   data = iris,
#'   regression_formula = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
#'   learners = learners[1:3])
#'
#' sl_model_predictions <- sl_model(iris)
#' fit_individual_learners <- lapply(learners[1:3], function(learner) { learner(data = iris, regression_formula = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width) } )
#' individual_learners_rmse <- lapply(fit_individual_learners, function(fit_learner) { rmse(fit_learner(iris) - iris$Sepal.Length) })
#'
#' print(paste0("super-learner rmse: ", rmse(sl_model_predictions - iris$Sepal.Length)))
#' individual_learners_rmse
#' }
#'
#' @export
super_learner <- function(
  data,
  learners,
  regression_formulas,
  y_variable,
  n_folds = 5,
  determine_superlearner_weights = determine_superlearner_weights_nnls,
  continuous_or_discrete = 'continuous',
  extra_learner_args = NULL,
  verbose_output = FALSE) {

  # throw an error if the learners are not a named list
  if (! is.list(learners) | length(unique(names(learners))) != length(learners)) {
    stop("The learners passed to lmpti::super_learner must have (unique) names.")
  }

  # add a folds column
  data <- make_folds(data, n_folds = n_folds)

  # split into an n_folds length list of training_data and validation_data:
  # training_data contains n_folds-1 folds of data, validation_data contains 1 fold
  training_data <- lapply(
    1:n_folds, function(i) {data |> dplyr::filter(sl_fold != i) |> dplyr::select(-sl_fold)})
  validation_data <- lapply(
    1:n_folds, function(i) {data |> dplyr::filter(sl_fold == i) |> dplyr::select(-sl_fold)})
  data$sl_fold <- NULL # remove sl_fold column now

  # make a tibble/dataframe to hold the trained learners:
  # one for each combination of a specific fold and a specific model
  trained_learners <- tibble::tibble(
    split = rep(1:n_folds, length(learners)),
    learner_name = rep(names(learners), each = n_folds))

  # handle vectorized regression_formulas argument
  #
  # if the regression_formulas is just a single formula, then we repeat it
  # in a vector length(learners) times to make it simple to just pass the ith
  # learner regression_formula[[i]].
  #
  if (inherits(regression_formulas, 'formula')) {
    regression_formulas <- rep(c(regression_formulas), length(learners)) # repeat the regression formula
  } else if (! (is.vector(regression_formulas) &&
                length(regression_formulas) == length(learners) &&
                all(sapply(regression_formulas, class) == 'formula'))) {
    stop("The regression_formula must either be a single formula or a vector
of formulas of the same length as the number of learners specified.")
  }

  # for each i in 1:n_folds and each model, train the model
  trained_learners$learned_predictor <- lapply(
    1:nrow(trained_learners), function(i) {
      # calculate which learner has the name for this row and use
      # the appropriate regression formula as well as the right
      # extra_learner_args
      learner_index <- which(names(learners) == trained_learners$learner_name[[i]])[[1]]

      # train the learner — the returned output is the prediction function from
      # the trained learner
      do.call(
        what = learners[[trained_learners[[i,'learner_name']]]],
        args = c(list(
          data = training_data[[trained_learners[[i,'split']]]],
          regression_formula = regression_formulas[[
            learner_index
          ]]),
          extra_learner_args[[learner_index]]
        )
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
  y_variables <- sapply(regression_formulas, function(f) as.character(f)[[2]])
  if (missing(y_variable) & length(unique(y_variables)) == 1) {
    y_variable <- unique(y_variables)
  } else if (missing(y_variable) & length(unique(y_variables)) > 1) {
    stop("Cannot infer the y-variable from the formulas passed.
Please pass yvar = ... to lmtp::super_learner.")
  }

  if (! y_variable %in% colnames(data)) {
    stop("The left-hand-side of the regression formula given must appear as a column in the data passed.")
  }

  # if the y_variable matches with any of the learners, we have problems —
  # the output second_stage_SL_dataset wouldn't be interpretable.
  if (y_variable %in% names(learners)) {
    stop("The outcome and names of all of the learners must be distinct, because the output
from super_learner is a data.frame with columns including the outcome variable and each of
the learners.")
  }

  # insert the validation Y data in another column next to the predictions
  second_stage_SL_dataset[[y_variable]] <- lapply(1:nrow(second_stage_SL_dataset), function(i) {
    validation_data[[second_stage_SL_dataset[[i, 'split']]]][[y_variable]]
  })

  # unnest all of the data (each cell prior to this contained a vector of either
  # predictions or the validation data)
  second_stage_SL_dataset <- tidyr::unnest(second_stage_SL_dataset, cols = colnames(second_stage_SL_dataset))

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
  learner_weights <- determine_superlearner_weights(second_stage_SL_dataset, y_variable)

  # adjust weights according to if using continuous or discrete super-learner
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
      do.call(
        what = learners[[i]],
        args = c(list(
          data = data, regression_formula = regression_formulas[[i]]
          ),
          extra_learner_args[[i]]
        )
      )
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

  if (verbose_output) {
    return(
      list(
        sl_predictor = predict_from_super_learned_model,
        holdout_predictions = second_stage_SL_dataset
        )
      )
  } else {
    return(predict_from_super_learned_model)
  }
}


#' Determine SuperLearner Weights with Nonnegative Least Squares
#'
#' This function accepts a dataframe that is structured to have
#' one column `Y` and other columns with unique names corresponding to
#' different model predictions for `Y`, and it will use nonnegative
#' least squares to determine the weights to use for a SuperLearner.
#'
#' @export
determine_superlearner_weights_nnls <- function(data, y_variable) {
  # use nonlinear least squares to produce a weighting scheme
  index_of_yvar <- which(colnames(data) == y_variable)[[1]]
  nnls_output <- nnls::nnls(
    A = as.matrix(data[,-index_of_yvar]),
    b = data[[y_variable]])

  return(nnls_output$x)
}

#' Add a new column called folds to the data passed
#'
#' A helper function for inserting a column called `sl_fold` into the
#' data with integers randomly sampled with replacement from `1:n_folds`.
#'
make_folds <- function(data, n_folds = 5) {
  if ('sl_fold' %in% colnames(data)) {
    stop("The data passed to make_folds already has a sl_fold column")
  }
  data$sl_fold <- sample.int(n = n_folds, size = nrow(data), replace = TRUE)
  return(data)
}
