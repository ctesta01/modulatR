# LMTP Data Structure ----------------------------------------------------------

#' Data structure for longitudinal LMTP analyses
#'
#' @description
#' `LMTP_Data_Struct` is a light-weight R6 container for longitudinal data
#' organized as baseline covariates `W`, time-varying covariates `L_t`,
#' treatments `A_t`, and a final outcome `Y`. It provides convenient methods
#' to extract per-time histories `H_t`, treatments, and outcomes needed by
#' the LMTP estimators.
#'
#' @export
LMTP_Data_Struct <- R6::R6Class(
  "LMTP_Data_Struct",
  public = list(
    #' @field data Data frame containing `id`, `W`, `L_t`, `A_t`, and `Y`.
    data = NULL,
    #' @field id_col Name of the id column.
    id_col = NULL,
    #' @field n_timesteps Number of treatment times (tau).
    n_timesteps = NULL,
    #' @field A_cols Character vector of treatment column names.
    A_cols = NULL,
    #' @field L_cols List of length `tau` of time-varying covariate names.
    L_cols = NULL,
    #' @field W_cols Optional character vector of baseline covariate names.
    W_cols = NULL,
    #' @field Y_col Outcome column name.
    Y_col = NULL,

    #' @description Construct a new LMTP data structure.
    #' @param data Data frame with all variables.
    #' @param id_col Name of id column.
    #' @param n_timesteps Integer number of treatment times.
    #' @param A_cols Character vector of treatment column names.
    #' @param L_cols List of character vectors for time-varying covariates.
    #' @param W_cols Optional character vector of baseline covariate names.
    #' @param Y_col Outcome column name.
    initialize = function(data, id_col, n_timesteps, A_cols, L_cols,
                          W_cols = NULL, Y_col) {
      self$data <- data
      self$id_col <- id_col
      self$n_timesteps <- as.integer(n_timesteps)
      self$A_cols <- A_cols
      self$L_cols <- L_cols
      self$W_cols <- W_cols
      self$Y_col <- Y_col
      private$validate()
    },

    #' @description Number of time points.
    #' @return Integer `tau`.
    tau = function() self$n_timesteps,

    #' @description Get the history `H_t = (W, L_1..L_t, A_1..A_{t-1})`.
    #' @param t Time index (1-based).
    #' @return Data frame of covariates in the history.
    H = function(t) {
      stopifnot(t >= 1, t <= self$tau())
      A_hist <- if (t > 1) self$A_cols[seq_len(t - 1)] else character(0)
      L_hist <- unlist(self$L_cols[seq_len(t)], use.names = FALSE)
      cols <- c(A_hist, L_hist, if (!is.null(self$W_cols)) self$W_cols)
      self$data[, cols, drop = FALSE]
    },

    #' @description Get `(A_t, H_t)` as a single data frame.
    #' @param t Time index (1-based).
    #' @return Data frame with `H_t` and the treatment column `A_t`.
    AH = function(t) {
      df <- self$H(t)
      df[[self$A_cols[[t]]]] <- self$data[[self$A_cols[[t]]]]
      df
    },

    #' @description Get the treatment vector at time `t`.
    #' @param t Time index.
    #' @return Vector `A_t`.
    A = function(t) self$data[[self$A_cols[[t]]]],

    #' @description Get the time-varying covariates at time `t`.
    #' @param t Time index.
    #' @return Data frame of covariates `L_t`.
    L = function(t) self$data[, self$L_cols[[t]], drop = FALSE],

    #' @description Get baseline covariates `W`.
    #' @return Data frame of baseline covariates or `NULL`.
    W = function() if (is.null(self$W_cols)) NULL else self$data[, self$W_cols, drop = FALSE],

    #' @description Get the outcome vector `Y`.
    #' @return Outcome vector.
    Y = function() self$data[[self$Y_col]],

    #' @description Convenience: return a list with `A`, `H`, and `H_next` at time `t`.
    #' @param t Time index.
    #' @return List with elements `A`, `H`, `AH`, and `H_next`.
    view_t = function(t) {
      list(
        A = self$A(t),
        H = self$H(t),
        AH = self$AH(t),
        H_next = if (t < self$tau()) self$H(t + 1) else NULL
      )
    }
  ),
  private = list(
    validate = function() {
      stopifnot(self$id_col %in% colnames(self$data))
      stopifnot(length(self$A_cols) == self$n_timesteps)
      stopifnot(is.list(self$L_cols), length(self$L_cols) == self$n_timesteps)
      stopifnot(all(vapply(self$L_cols, function(v) is.character(v) && length(v) >= 1, TRUE)))
      needed <- c(unlist(self$L_cols), self$A_cols,
                  if (!is.null(self$W_cols)) self$W_cols,
                  self$Y_col)
      if (!all(needed %in% colnames(self$data))) {
        miss <- setdiff(needed, colnames(self$data))
        stop("LMTP_Data_Struct: missing columns in data: ", paste(miss, collapse = ", "))
      }
    }
  )
)
