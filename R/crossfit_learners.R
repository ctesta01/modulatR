
#' Cross-fit an LMTP learner
#'
#' @description
#' Wraps any LMTP learner in an origami V-fold cross-fitting layer.
#'
#' The input `learner` must take one `LMTPData` object and return a prediction
#' function. The wrapped learner has the same external contract:
#'
#'   crossfit_lmtp_learner(learner)(data) -> function(newdata)
#'
#' For rows whose ids are present in the training object, predictions are made
#' with the fold-specific model that did not train on that row. For new rows not
#' present in the original training object, predictions are averaged over the
#' fold-specific models.
#'
#' @export
# Cross-fitting wrapper ----------------------------------------------------

crossfit_lmtp_learner <- function(learner,
                                  V = 5,
                                  folds = NULL,
                                  use_future = FALSE,
                                  fold_fun = origami::folds_vfold,
                                  ...) {
  if (!is.function(learner)) {
    stop("`learner` must be a function.")
  }

  force(learner)
  force(V)
  force(folds)
  force(use_future)
  force(fold_fun)

  function(data) {
    if (!inherits(data, "LMTPData")) {
      stop("`data` must inherit from `LMTPData`.")
    }

    folds_local <- folds %||% origami::make_folds(
      n = data$n,
      fold_fun = fold_fun,
      V = V,
      ...
    )

    get_fold_indices <- function(fold, which = c("training", "validation")) {
      which <- match.arg(which)

      candidates <- switch(
        which,
        training = c("training", "training_set", "train", "train_set"),
        validation = c("validation", "validation_set", "valid", "valid_set")
      )

      out <- NULL

      for (nm in candidates) {
        val <- tryCatch(fold[[nm]], error = function(e) NULL)
        if (!is.null(val)) {
          out <- val
          break
        }
      }

      if (is.null(out)) {
        stop(
          "Could not find a `", which, "` index field in an origami fold. ",
          "Available fields are: ",
          paste(names(fold), collapse = ", ")
        )
      }

      if (is.logical(out)) {
        out <- which(out)
      }

      as.integer(out)
    }

    train_one_fold <- function(j) {
      fold <- folds_local[[j]]

      train_idx <- get_fold_indices(fold, "training")
      valid_idx <- get_fold_indices(fold, "validation")

      data_train <- .subset_ds(data, train_idx)
      predict_train <- learner(data_train)

      if (!is.function(predict_train)) {
        stop("The base learner must return a prediction function.")
      }

      list(
        fold_index = j,
        training = train_idx,
        validation = valid_idx,
        predict = predict_train
      )
    }

    fold_results <- if (isTRUE(use_future)) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        stop("`future.apply` is required when `use_future = TRUE`.")
      }
      future.apply::future_lapply(seq_along(folds_local), train_one_fold)
    } else {
      lapply(seq_along(folds_local), train_one_fold)
    }

    id_col <- data$id_col
    id_full <- data$df[[id_col]]

    validation_lookup <- rep(NA_integer_, data$n)

    for (j in seq_along(fold_results)) {
      validation_lookup[fold_results[[j]]$validation] <- j
    }

    if (anyNA(validation_lookup)) {
      bad <- which(is.na(validation_lookup))
      stop(
        "Some rows were not assigned to a validation fold. ",
        "First missing rows: ",
        paste(utils::head(bad, 10), collapse = ", ")
      )
    }

    function(newdata) {
      if (!inherits(newdata, "LMTPData")) {
        stop("`newdata` must inherit from `LMTPData`.")
      }

      new_ids <- newdata$df[[id_col]]
      idx_in_original <- match(new_ids, id_full)

      preds <- rep(NA_real_, newdata$n)

      in_original <- which(!is.na(idx_in_original))

      if (any(! in_original)) {
        warning("newdata contained observations that were not used in the cross-fitting of nuisance learners")
        warning("rows: ", paste0(which( ! in_original ), collapse = ', '))
      }

      if (length(in_original) > 0L) {
        fold_for_row <- validation_lookup[idx_in_original[in_original]]

        for (j in seq_along(fold_results)) {
          row_pos <- in_original[fold_for_row == j]

          if (length(row_pos) == 0L) next

          preds[row_pos] <- as.numeric(
            fold_results[[j]]$predict(.subset_ds(newdata, row_pos))
          )
        }
      }

      new_rows <- which(is.na(idx_in_original))

      if (length(new_rows) > 0L) {
        pred_mat <- vapply(
          fold_results,
          function(fr) {
            as.numeric(fr$predict(.subset_ds(newdata, new_rows)))
          },
          numeric(length(new_rows))
        )

        preds[new_rows] <- rowMeans(pred_mat)
      }

      if (anyNA(preds)) {
        stop("Cross-fitted prediction failed for at least one row.")
      }

      preds
    }
  }
}
