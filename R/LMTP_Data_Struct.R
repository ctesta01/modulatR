#' Longitudinal data object for LMTP analyses
#'
#' @description
#' `LMTPData` is an R6 class that provides structure to longitudinal causal data
#' organized as baseline covariates `W`,
#' time-varying covariates `L_t`, treatments `A_t`,
#' and a final outcome `Y`.
#'
#' It provides convenient methods to extract per-time histories `H_t`,
#' treatments, outcomes, and related longitudinal views used throughout
#' `{modulatR}`.
#'
#' The object is intentionally structural rather than algorithmic: it stores
#' the node layout and exposes standardized views of the observed data, but it
#' does not itself perform nuisance estimation, imputation, or targeting.
#'
#' @export
#' @examples
#' ## Single timepoint example (tau = 1)
#' n <- 100
#' df <- tibble::tibble(
#'   `L[t1,1]` = rnorm(n), `L[t1,2]` = rbinom(n, 1, 0.5),
#'   A = rbinom(n, 1, plogis(0.4 + `L[t1,1]` + `L[t1,2]`)),
#'   Y = rnorm(n, mean = `L[t1,2]` + 2.3*A, sd = .1))
#' ds <- LMTPData$new(
#'   data = df,
#'   L_cols = list(c('L[t1,2]', 'L[t1,2]')),
#'   A_cols = 'A',
#'   Y_col = 'Y')
#' ds
#'
#' ## Two-timepoint example (tau = 2)
#' df <- data.frame(
#'   id = 1:n,
#'   W1 = rnorm(n),
#'   L1 = rnorm(n),
#'   A1 = rbinom(n, 1, 0.5),
#'   L2 = rnorm(n),
#'   A2 = rbinom(n, 1, 0.5),
#'   Y  = rnorm(n)
#' )
#'
#' ds <- LMTPData$new(
#'   data = df,
#'   id_col = "id",
#'   A_cols = c("A1", "A2"),
#'   L_cols = list(
#'     "L1",
#'     "L2"
#'   ),
#'   W_cols = "W1",
#'   Y_col = "Y"
#' )
#'
#' ds
#'
#' # Histories
#' head(ds$H(1))  # (W, L1)
#' head(ds$H(2))  # (W, L1, L2, A1)
#'
#' # Treatment at time 2
#' head(ds$A(2))
#' head(ds$AH(2))
LMTPData <- R6::R6Class(
  "LMTPData",
  lock_objects = FALSE,
  public = list(
    #' @field id_col Name of the unit identifier column.
    id_col = NULL,

    #' @field A_cols Character vector of treatment column names.
    A_cols = NULL,

    #' @field L_cols List of character vectors giving time-varying covariates
    #'   at each time point.
    L_cols = NULL,

    #' @field W_cols Optional character vector of baseline covariate names.
    W_cols = NULL,

    #' @field Y_col Outcome column name.
    Y_col = NULL,

    #' @field metadata Optional list of additional metadata.
    metadata = NULL,

    #' @description Create a new LMTPData object.
    #' @param data A data.frame containing all variables.
    #' @param id_col Name of the unit identifier column.
    #' @param W_cols Optional character vector of baseline covariate names.
    #' @param L_cols List of character vectors for time-varying covariates.
    #' @param A_cols Character vector of treatment column names, one per time.
    #' @param Y_col Outcome column name.
    #' @param metadata Optional list of metadata.
    initialize = function(data,
                          id_col = NULL,
                          W_cols = NULL,
                          L_cols = NULL,
                          A_cols,
                          Y_col,
                          metadata = list()) {
      private$.data <- as.data.frame(data)

      if (is.null(id_col)) {
        auto_id <- "..row_id"
        if (auto_id %in% colnames(private$.data)) {
          stop(
            "No `id_col` was supplied, but the reserved auto-generated id column ",
            "`", auto_id, "` already exists in `data`. Please supply `id_col` explicitly."
          )
        }
        private$.data[[auto_id]] <- seq_len(nrow(private$.data))
        message("`id_col` was NULL; using `", auto_id, "` = 1:nrow(data).")
        id_col <- auto_id
      }

      self$id_col <- id_col
      self$A_cols <- A_cols
      self$L_cols <- L_cols
      self$W_cols <- W_cols
      self$Y_col <- Y_col
      self$metadata <- metadata

      private$validate()
      invisible(self)
    },

    #' @description Number of treatment times.
    #' @return Integer number of treatment times.
    tau = function() {
      length(self$A_cols)
    },

    #' @description Return the treatment vector `A_t`.
    #' @param t Time index.
    #' @return A vector.
    A = function(t) {
      private$check_t(t)
      private$.data[[self$A_cols[[t]]]]
    },

    #' @description Return the time-varying covariates `L_t`.
    #' @param t Time index.
    #' @return A data.frame.
    L = function(t) {
      private$check_t(t)
      private$.data[, self$L_cols[[t]], drop = FALSE]
    },

    #' @description Return the baseline covariates `W`.
    #' @return A data.frame or `NULL`.
    W = function() {
      if (is.null(self$W_cols) || length(self$W_cols) == 0L) {
        return(NULL)
      }
      private$.data[, self$W_cols, drop = FALSE]
    },

    #' @description Return the outcome vector `Y`.
    #' @return A vector.
    Y = function() {
      private$.data[[self$Y_col]]
    },

    #' @description Return the history `H_t = (W, \bar{L}_t, \bar{A}_{t-1})`.
    #' @param t Time index.
    #' @return A data.frame.
    H = function(t) {
      cols <- private$ordered_history_cols(t, include_A_t = FALSE, include_Y = FALSE)
      private$.data[, cols, drop = FALSE]
      },

    #' @description Return `(H_t, A_t)` as one data frame.
    #' @param t Time index.
    #' @return A data.frame.
    AH = function(t) {
      cols <- private$ordered_history_cols(t, include_A_t = TRUE, include_Y = FALSE)
      private$.data[, cols, drop = FALSE]
    },

    #' @description Return the next history `H_{t+1}` when it exists.
    #' @param t Time index.
    #' @return A data.frame or `NULL`.
    H_next = function(t) {
      private$check_t(t)
      if (t >= self$tau()) {
        return(NULL)
      }
      self$H(t + 1L)
    },

    #' @description Return all node names used by the object.
    #' @return A named list.
    nodes = function() {
      list(
        id = self$id_col,
        W = self$W_cols,
        L = self$L_cols,
        A = self$A_cols,
        Y = self$Y_col
      )
    },

    #' @description Print a compact summary of the data object.
    #' @return The object invisibly.
    print = function(...) {
      cat("LMTPData\n")
      cat("  n: ", self$n, "\n", sep = "")
      cat("  tau: ", self$tau(), "\n", sep = "")
      if (! is.null(self$id_col)) { cat("  id: ", self$id_col, "\n", sep = "") }

      if (!is.null(self$W_cols) && length(self$W_cols) > 0L) {
        cat("  W: ", paste(self$W_cols, collapse = ", "), "\n", sep = "")
      } else {
        cat("  W: <none>\n")
      }

      cat("  A: ", paste(self$A_cols, collapse = ", "), "\n", sep = "")
      cat("  Y: ", self$Y_col, "\n", sep = "")

      cat("  L by time:\n")
      for (t in seq_len(self$tau())) {
        cat("    t = ", t, ": ", paste(self$L_cols[[t]], collapse = ", "), "\n", sep = "")
      }

      invisible(self)
    }
  ),

  active = list(
    #' @description Return the underlying data frame.
    #' @return A data.frame.
    df = function(value) {
      if (!missing(value)) stop("`df` is read-only.")
      private$.data
    },

    #' @field n Number of rows in the data.
    n = function(value) {
      if (!missing(value)) stop("`n` is read-only.")
      nrow(private$.data)
    },

    #' @field p Number of columns in the underlying data.
    p = function(value) {
      if (!missing(value)) stop("`p` is read-only.")
      ncol(private$.data)
    },

    #' @field colnames Column names of the underlying data.
    colnames = function(value) {
      if (!missing(value)) stop("`colnames` is read-only.")
      colnames(private$.data)
    }
  ),

  private = list(
    .data = NULL,

    check_t = function(t) {
      if (!is.numeric(t) || length(t) != 1L || is.na(t)) {
        stop("`t` must be a single non-missing numeric value.")
      }

      t <- as.integer(t)

      if (t < 1L || t > length(self$A_cols)) {
        stop("`t` must be between 1 and tau = ", length(self$A_cols), ".")
      }

      invisible(t)
    },

    validate = function() {
      if (!is.data.frame(private$.data)) {
        stop("`data` must be a data.frame.")
      }

      if (!is.null(self$id_col)) {
        if (!is.character(self$id_col) || length(self$id_col) != 1L) {
          stop("`id_col` must be a single character string.")
        }
        if (!self$id_col %in% colnames(private$.data)) {
          stop("`id_col` not found in data: ", self$id_col)
        }
      }


      if (!is.character(self$A_cols) || length(self$A_cols) < 1L) {
        stop("`A_cols` must be a non-empty character vector.")
      }

      if (!is.list(self$L_cols) || length(self$L_cols) != length(self$A_cols)) {
        stop("`L_cols` must be a list of the same length as `A_cols`.")
      }

      ok_L <- vapply(
        self$L_cols,
        function(x) is.character(x) && length(x) >= 1L,
        logical(1)
      )

      if (!all(ok_L)) {
        stop("Each element of `L_cols` must be a non-empty character vector.")
      }

      if (!is.null(self$W_cols) && !is.character(self$W_cols)) {
        stop("`W_cols` must be NULL or a character vector.")
      }

      if (!is.character(self$Y_col) || length(self$Y_col) != 1L) {
        stop("`Y_col` must be a single character string.")
      }

      needed <- c(
        self$id_col,
        self$A_cols,
        unlist(self$L_cols, use.names = FALSE),
        if (!is.null(self$W_cols)) self$W_cols,
        self$Y_col
      )

      missing_cols <- setdiff(unique(needed), colnames(private$.data))
      if (length(missing_cols) > 0L) {
        stop(
          "LMTPData is missing required columns: ",
          paste(missing_cols, collapse = ", ")
        )
      }
    },

    ordered_history_cols = function(t, include_A_t = FALSE, include_Y = FALSE) {
      private$check_t(t)

      cols <- character(0)

      if (!is.null(self$id_col)) {
        cols <- c(cols, self$id_col)
      }

      if (!is.null(self$W_cols) && length(self$W_cols) > 0L) {
        cols <- c(cols, self$W_cols)
      }

      for (s in seq_len(t)) {
        cols <- c(cols, self$L_cols[[s]])

        if (s < t || include_A_t) {
          cols <- c(cols, self$A_cols[[s]])
        }
      }

      if (include_Y) {
        cols <- c(cols, self$Y_col)
      }

      cols
    }
  )
)
