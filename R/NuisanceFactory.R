
#' LMTP nuisance factory
#'
#' @description
#' `LMTPNuisanceFactory` coordinates nuisance estimation for LMTP estimators.
#'
#' Users provide time-indexed learner functions. Each learner takes a single
#' `LMTPData` object and returns a prediction function. The factory uses these
#' prediction functions immediately, stores prediction vectors, and does not keep
#' prediction functions as trained artifacts.
#'
#' Public learner contracts:
#'
#' m learner:
#'   m_learner(data_t) -> function(newdata_t) numeric(nrow(newdata_t$df))
#'
#' g learner:
#'   g_learner(data_t) -> function(newdata_t) numeric(nrow(newdata_t$df))
#'
#  lambda learner:
#'   lambda_learner(data_t_augmented) -> function(newdata_t) numeric(nrow(newdata_t$df))
#'
#' A m learner predicts the pseudo-outcome stored in
#' `data_t$metadata$pseudo_outcome_col`.
#'
#' A g learner predicts the observed treatment density or mass
#' `g_t(A_t | H_t)`.
#'
#' A lambda learner predicts `P(lambda = 1 | A_t, H_t)` in the augmented
#' observed-versus-policy-shifted classification problem where `lambda = 1` is the
#' indicator for when `A_t` is the intervention value. The factory converts
#' this to the density ratio `r_t = lambda_t / (1 - lambda_t)`.
#'
#' Exactly one of `g_learners` and `lambda_learners` must be supplied.
#'
#' @export
# LMTP nuisance factory ----------------------------------------------------

# LMTP nuisance factory ----------------------------------------------------

LMTPNuisanceFactory <- R6::R6Class(
  "LMTPNuisanceFactory",
  public = list(
    tau = NULL,
    policy_seq = NULL,

    m_learners = NULL,
    g_learners = NULL,
    lambda_learners = NULL,

    truncate_density = NULL,
    truncate_ratio = NULL,
    clip_lambda_probability = NULL,

    m_obs_preds = NULL,
    m_d_preds = NULL,

    g_obs_preds = NULL,
    gd_obs_preds = NULL,
    g_d_input_preds = NULL,
    gd_d_input_preds = NULL,

    lambda_obs_preds = NULL,
    lambda_d_preds = NULL,

    r_obs_preds = NULL,
    r_d_preds = NULL,
    omega_preds = NULL,

    initialize = function(policy_seq,
                          m_learners,
                          g_learners = NULL,
                          lambda_learners = NULL,
                          truncate_density = 1e-12,
                          truncate_ratio = c(1e-6, Inf),
                          clip_lambda_probability = 1e-6) {
      if (!inherits(policy_seq, "LMTPPolicySequence")) {
        stop("`policy_seq` must inherit from `LMTPPolicySequence`.")
      }

      self$tau <- policy_seq$tau()
      self$policy_seq <- policy_seq

      self$m_learners <- private$as_time_list(
        m_learners,
        self$tau,
        "m_learners"
      )

      has_g <- !is.null(g_learners)
      has_lambda <- !is.null(lambda_learners)

      if (has_g == has_lambda) {
        stop("Supply exactly one of `g_learners` or `lambda_learners`.")
      }

      self$g_learners <- if (has_g) {
        private$as_time_list(g_learners, self$tau, "g_learners")
      } else {
        NULL
      }

      self$lambda_learners <- if (has_lambda) {
        private$as_time_list(lambda_learners, self$tau, "lambda_learners")
      } else {
        NULL
      }

      if (!is.numeric(truncate_density) ||
          length(truncate_density) != 1L ||
          truncate_density <= 0) {
        stop("`truncate_density` must be a positive scalar.")
      }

      if (!is.numeric(truncate_ratio) ||
          length(truncate_ratio) != 2L ||
          truncate_ratio[1] <= 0 ||
          truncate_ratio[2] < truncate_ratio[1]) {
        stop("`truncate_ratio` must be c(lower, upper), with 0 < lower <= upper.")
      }

      if (!is.numeric(clip_lambda_probability) ||
          length(clip_lambda_probability) != 1L ||
          clip_lambda_probability <= 0 ||
          clip_lambda_probability >= 0.5) {
        stop("`clip_lambda_probability` must be a scalar in (0, 0.5).")
      }

      self$truncate_density <- truncate_density
      self$truncate_ratio <- truncate_ratio
      self$clip_lambda_probability <- clip_lambda_probability

      invisible(self)
    },

    validate_ds = function(ds) {
      if (!inherits(ds, "LMTPData")) {
        stop("`ds` must inherit from `LMTPData`.")
      }

      if (ds$tau() != self$tau) {
        stop(
          "`ds$tau()` is ", ds$tau(),
          ", but the nuisance factory has tau = ", self$tau, "."
        )
      }

      self$policy_seq$validate_against_data(ds)
      invisible(TRUE)
    },

    train_m_t = function(ds,
                         t,
                         pseudo_outcome,
                         pseudo_outcome_col = "..pseudo_outcome") {
      self$validate_ds(ds)
      private$check_t(t)

      ds_t <- private$prepare_time_data(ds, t)
      ds_t <- .with_pseudo_outcome(
        data = ds_t,
        pseudo_outcome = pseudo_outcome,
        pseudo_outcome_col = pseudo_outcome_col
      )

      m_predict <- self$m_learners[[t]](ds_t)

      if (!is.function(m_predict)) {
        stop("`m_learners[[", t, "]]` must return a prediction function.")
      }

      ds_t_d <- private$with_modified_treatment(ds_t, t)

      m_obs <- private$as_numeric_prediction(m_predict(ds_t), ds_t$n, "m_obs")
      m_d <- private$as_numeric_prediction(m_predict(ds_t_d), ds_t$n, "m_d")

      if (is.null(self$m_obs_preds)) self$m_obs_preds <- vector("list", self$tau)
      if (is.null(self$m_d_preds)) self$m_d_preds <- vector("list", self$tau)

      self$m_obs_preds[[t]] <- m_obs
      self$m_d_preds[[t]] <- m_d

      list(
        m_obs = m_obs,
        m_d = m_d
      )
    },

    train_m_models = function(ds,
                              pseudo_outcome_col = "..pseudo_outcome",
                              quiet = FALSE) {
      self$validate_ds(ds)

      if (!isTRUE(quiet)) {
        message(
          "`train_m_models()` fits the sequential outcome-regression recursion. ",
          "For the targeted estimator, use `run_tmle_for_LMTP()`."
        )
      }

      self$m_obs_preds <- vector("list", self$tau)
      self$m_d_preds <- vector("list", self$tau)

      pseudo_outcome <- ds$Y()

      for (t in rev(seq_len(self$tau))) {
        m_t <- self$train_m_t(
          ds = ds,
          t = t,
          pseudo_outcome = pseudo_outcome,
          pseudo_outcome_col = pseudo_outcome_col
        )

        if (t > 1L) {
          pseudo_outcome <- m_t$m_d
        }
      }

      invisible(self)
    },

    train_ratios = function(ds) {
      self$validate_ds(ds)

      self$r_obs_preds <- vector("list", self$tau)
      self$r_d_preds <- vector("list", self$tau)
      self$omega_preds <- vector("list", self$tau)

      if (!is.null(self$g_learners)) {
        self$g_obs_preds <- vector("list", self$tau)
        self$gd_obs_preds <- vector("list", self$tau)
        self$g_d_input_preds <- vector("list", self$tau)
        self$gd_d_input_preds <- vector("list", self$tau)
      }

      if (!is.null(self$lambda_learners)) {
        self$lambda_obs_preds <- vector("list", self$tau)
        self$lambda_d_preds <- vector("list", self$tau)
      }

      omega_t <- rep(1, ds$n)

      for (t in seq_len(self$tau)) {
        ratio_t <- if (!is.null(self$g_learners)) {
          self$train_g_t(ds, t)
        } else {
          self$train_lambda_t(ds, t)
        }

        self$r_obs_preds[[t]] <- .truncate_interval(
          ratio_t$r_obs,
          truncate = self$truncate_ratio
        )

        self$r_d_preds[[t]] <- .truncate_interval(
          ratio_t$r_d,
          truncate = self$truncate_ratio
        )

        omega_t <- omega_t * self$r_obs_preds[[t]]
        self$omega_preds[[t]] <- omega_t
      }

      invisible(self)
    },

    train = function(ds) {
      self$train_ratios(ds)
    },

    train_g_t = function(ds, t) {
      self$validate_ds(ds)
      private$check_t(t)

      if (is.null(self$g_learners)) {
        stop("`g_learners` is NULL.")
      }

      ds_t <- private$prepare_time_data(ds, t)
      ds_t_d <- private$with_modified_treatment(ds_t, t)

      g_predict <- self$g_learners[[t]](ds_t)

      if (!is.function(g_predict)) {
        stop("`g_learners[[", t, "]]` must return a prediction function.")
      }

      A_t <- ds_t$A(t)
      H_t <- ds_t$H(t)
      A_t_d <- ds_t_d$A(t)

      g_obs <- .truncate_positive(
        private$as_numeric_prediction(g_predict(ds_t), ds_t$n, "g_obs"),
        truncate_density = self$truncate_density
      )

      g_d_input <- .truncate_positive(
        private$as_numeric_prediction(g_predict(ds_t_d), ds_t$n, "g_d_input"),
        truncate_density = self$truncate_density
      )

      density_fun <- function(A_vec, H_df) {
        newdata <- private$make_prediction_data(
          template = ds_t,
          t = t,
          A_vec = A_vec,
          H_df = H_df
        )

        private$as_numeric_prediction(
          g_predict(newdata),
          length(A_vec),
          "density_fun"
        )
      }

      gd_obs <- .truncate_positive(
        self$policy_seq$gd_t(
          t = t,
          A_vec = A_t,
          H_df = H_t,
          density_fun = density_fun
        ),
        truncate_density = self$truncate_density
      )

      gd_d_input <- .truncate_positive(
        self$policy_seq$gd_t(
          t = t,
          A_vec = A_t_d,
          H_df = H_t,
          density_fun = density_fun
        ),
        truncate_density = self$truncate_density
      )

      self$g_obs_preds[[t]] <- g_obs
      self$g_d_input_preds[[t]] <- g_d_input
      self$gd_obs_preds[[t]] <- gd_obs
      self$gd_d_input_preds[[t]] <- gd_d_input

      list(
        r_obs = gd_obs / g_obs,
        r_d = gd_d_input / g_d_input
      )
    },

    train_lambda_t = function(ds, t) {
      self$validate_ds(ds)
      private$check_t(t)

      if (is.null(self$lambda_learners)) {
        stop("`lambda_learners` is NULL.")
      }

      ds_t <- private$prepare_time_data(ds, t)
      ds_t_d <- private$with_modified_treatment(ds_t, t)

      ds_aug <- private$augment_for_ratio_classification(ds_t, t)

      lambda_predict <- self$lambda_learners[[t]](ds_aug)

      if (!is.function(lambda_predict)) {
        stop("`lambda_learners[[", t, "]]` must return a prediction function.")
      }

      lambda_obs <- .clip_probability(
        private$as_numeric_prediction(
          lambda_predict(ds_t),
          ds_t$n,
          "lambda_obs"
        ),
        clip_probability = self$clip_lambda_probability
      )

      lambda_d <- .clip_probability(
        private$as_numeric_prediction(
          lambda_predict(ds_t_d),
          ds_t$n,
          "lambda_d"
        ),
        clip_probability = self$clip_lambda_probability
      )

      self$lambda_obs_preds[[t]] <- lambda_obs
      self$lambda_d_preds[[t]] <- lambda_d

      list(
        r_obs = lambda_obs / (1 - lambda_obs),
        r_d = lambda_d / (1 - lambda_d)
      )
    },

    omega = function(t) {
      private$check_ratios_trained()
      private$check_t(t)
      self$omega_preds[[t]]
    },

    omega_prev = function(t, n = NULL) {
      private$check_ratios_trained()
      private$check_t(t)

      if (t == 1L) {
        if (is.null(n)) n <- length(self$omega_preds[[1L]])
        return(rep(1, n))
      }

      self$omega_preds[[t - 1L]]
    },

    print = function(...) {
      ratio_source <- if (!is.null(self$g_learners)) {
        "g_learners"
      } else {
        "lambda_learners"
      }

      cat("LMTPNuisanceFactory\n")
      cat("  tau: ", self$tau, "\n", sep = "")
      cat("  policy: ", self$policy_seq$name, "\n", sep = "")
      cat("  ratio source: ", ratio_source, "\n", sep = "")
      cat("  ratios trained: ", !is.null(self$r_obs_preds), "\n", sep = "")
      cat("  m regressions trained: ", !is.null(self$m_d_preds), "\n", sep = "")
      invisible(self)
    }
  ),

  private = list(
    as_time_list = function(x, tau, label) {
      if (is.function(x)) {
        return(rep(list(x), tau))
      }

      if (is.list(x) &&
          length(x) == tau &&
          all(vapply(x, is.function, logical(1)))) {
        return(x)
      }

      stop("`", label, "` must be a function or a list of functions of length tau.")
    },

    check_t = function(t) {
      if (!is.numeric(t) || length(t) != 1L || is.na(t)) {
        stop("`t` must be a single non-missing numeric value.")
      }

      t <- as.integer(t)

      if (t < 1L || t > self$tau) {
        stop("`t` must be between 1 and tau = ", self$tau, ".")
      }

      invisible(TRUE)
    },

    check_ratios_trained = function() {
      if (is.null(self$r_obs_preds) ||
          is.null(self$r_d_preds) ||
          is.null(self$omega_preds)) {
        stop("Call `$train_ratios(ds)` before requesting ratio or omega vectors.")
      }

      invisible(TRUE)
    },

    prepare_time_data = function(ds, t) {
      out <- .with_lmtp_time(ds, t)
      out
    },

    with_modified_treatment = function(ds_t, t) {
      A_t_d <- self$policy_seq$apply_t(t, ds_t$A(t), ds_t$H(t))

      private$replace_treatment_in_time_data(
        ds_t = ds_t,
        t = t,
        A_vec = A_t_d
      )
    },

    replace_treatment_in_time_data = function(ds_t, t, A_vec) {
      if (length(A_vec) != ds_t$n) {
        stop("`A_vec` must have length `ds_t$n`.")
      }

      df <- ds_t$df
      df[[ds_t$A_cols[[t]]]] <- A_vec

      LMTPData$new(
        data = df,
        id_col = ds_t$id_col,
        A_cols = ds_t$A_cols,
        L_cols = ds_t$L_cols,
        W_cols = ds_t$W_cols,
        Y_col = ds_t$Y_col,
        metadata = ds_t$metadata
      )
    },

    make_prediction_data = function(template, t, A_vec, H_df) {
      if (length(A_vec) != nrow(H_df)) {
        stop("`A_vec` must have length `nrow(H_df)`.")
      }

      df <- template$df[seq_len(length(A_vec)), , drop = FALSE]

      history_cols <- colnames(template$H(t))
      missing_history <- setdiff(history_cols, colnames(H_df))

      if (length(missing_history) > 0L) {
        stop(
          "`H_df` is missing required history columns: ",
          paste(missing_history, collapse = ", ")
        )
      }

      df[, history_cols] <- H_df[, history_cols, drop = FALSE]
      df[[template$A_cols[[t]]]] <- A_vec

      LMTPData$new(
        data = df,
        id_col = template$id_col,
        A_cols = template$A_cols,
        L_cols = template$L_cols,
        W_cols = template$W_cols,
        Y_col = template$Y_col,
        metadata = template$metadata
      )
    },

    augment_for_ratio_classification = function(ds_t,
                                                t,
                                                lambda_col = "..lambda") {
      private$check_t(t)

      A_t_d <- self$policy_seq$apply_t(
        t = t,
        A_vec = ds_t$A(t),
        H_df = ds_t$H(t)
      )

      df_obs <- ds_t$df
      df_obs[[lambda_col]] <- 0L

      df_shift <- ds_t$df
      df_shift[[ds_t$A_cols[[t]]]] <- A_t_d
      df_shift[[lambda_col]] <- 1L

      df_aug <- rbind(df_obs, df_shift)

      out <- LMTPData$new(
        data = df_aug,
        id_col = ds_t$id_col,
        A_cols = ds_t$A_cols,
        L_cols = ds_t$L_cols,
        W_cols = ds_t$W_cols,
        Y_col = ds_t$Y_col,
        metadata = ds_t$metadata
      )

      out$metadata$lambda_col <- lambda_col
      out$metadata$t <- t
      out
    },

    as_numeric_prediction = function(x, n, label) {
      x <- as.numeric(x)

      if (length(x) != n) {
        stop("`", label, "` must have length ", n, ".")
      }

      if (anyNA(x)) {
        stop("`", label, "` contains missing values.")
      }

      x
    }
  )
)
