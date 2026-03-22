#' LMTP fluctuation submodel
#'
#' @description
#' `LMTPFluctuationSubmodel` implements one-step TMLE fluctuations for
#' longitudinal modified treatment policy estimators.
#'
#' Unlike the generic fluctuation object in `{TargetedLearning}`, this class is
#' specialized for the LMTP recursion. In particular, it returns a wrapped
#' updated regression function `Q_t^*(A, H)` and supports clever covariates that
#' are either scalar- or vector-valued.
#'
#' The update is:
#'
#' * identity fluctuation:
#'   \deqn{Q_t^*(A,H) = Q_t(A,H) + H_t(A,H)^\top \epsilon}
#'
#' * logit fluctuation:
#'   \deqn{Q_t^*(A,H) =
#'     \mathrm{from01}\left[
#'       \mathrm{expit}\left(
#'         \mathrm{logit}(\mathrm{to01}(Q_t(A,H))) + H_t(A,H)^\top \epsilon
#'       \right)
#'     \right]}
#'
#' This class deliberately does not implement the Wei et al. Lagrangian-dual
#' stabilization; that should live in a separate fluctuation or targeting object.
#'
#' @export
LMTPFluctuationSubmodel <- R6::R6Class(
  "LMTPFluctuationSubmodel",
  public = list(
    family = NULL,
    bounds = c(NA_real_, NA_real_),
    clip = 1e-6,
    use_intercept = FALSE,
    to01 = NULL,
    from01 = NULL,

    #' @description Create a new LMTP fluctuation submodel.
    #' @param family Either `stats::gaussian()` for identity fluctuation or
    #'   `stats::binomial()` for logit fluctuation.
    #' @param bounds Optional lower/upper bounds for the outcome when using a
    #'   logit fluctuation. Defaults to `c(0, 1)`.
    #' @param clip Numeric truncation level for probabilities in logit updates.
    #' @param use_intercept Logical; whether to include an intercept in the
    #'   fluctuation regression.
    initialize = function(family = stats::gaussian(),
                          bounds = NULL,
                          clip = 1e-6,
                          use_intercept = FALSE) {
      is_logistic <- inherits(family, "family") &&
        family$family == "binomial" &&
        family$link == "logit"

      is_identity <- inherits(family, "family") &&
        family$family == "gaussian" &&
        family$link == "identity"

      if (!(is_identity || is_logistic)) {
        stop(
          "At this time, `LMTPFluctuationSubmodel` only supports ",
          "`gaussian(identity)` and `binomial(logit)` fluctuations."
        )
      }

      self$family <- family
      self$clip <- clip
      self$use_intercept <- use_intercept

      if (is_logistic) {
        if (is.null(bounds)) {
          bounds <- c(0, 1)
        }
        if (!is.numeric(bounds) || length(bounds) != 2L) {
          stop("`bounds` must be a numeric vector of length 2.")
        }
        bounds <- sort(bounds)
        self$bounds <- bounds

        self$to01 <- function(x) {
          (x - self$bounds[1]) / (self$bounds[2] - self$bounds[1])
        }
        self$from01 <- function(z) {
          self$bounds[1] + (self$bounds[2] - self$bounds[1]) * z
        }
      }

      invisible(self)
    },

    #' @description Fit a single TMLE fluctuation update.
    #'
    #' @param Q0_fun_t Function `(A_vec, H_df) -> numeric` returning the current
    #'   regression `Q_t(A,H)`.
    #' @param Qvec Numeric vector of current predictions at the observed data.
    #' @param target_vec Numeric vector of targets, e.g. `Y` or a next-stage
    #'   pseudo-outcome.
    #' @param H_obs Observed clever covariate values at the observed data.
    #'   May be a numeric vector, matrix, or data.frame.
    #' @param H_fun_t Function `(A_vec, H_df) -> H(A,H)` returning the clever
    #'   covariate at arbitrary evaluation points. This may return a vector,
    #'   matrix, or data.frame.
    #' @param t Optional time index used only in messages/errors.
    #'
    #' @return A list with components:
    #'   * `epsilon`: scalar or vector fluctuation coefficient(s),
    #'   * `intercept`: scalar intercept (0 unless `use_intercept = TRUE`),
    #'   * `fit`: fitted GLM object,
    #'   * `wrap`: updated regression function `(A_vec, H_df) -> Q_t^*(A,H)`.
    fit_update = function(Q0_fun_t,
                          Qvec,
                          target_vec,
                          H_obs,
                          H_fun_t,
                          t = NULL) {

      if (!is.function(Q0_fun_t)) {
        stop("`Q0_fun_t` must be a function.")
      }
      if (!is.function(H_fun_t)) {
        stop("`H_fun_t` must be a function.")
      }
      if (!is.numeric(Qvec) || !is.numeric(target_vec)) {
        stop("`Qvec` and `target_vec` must be numeric.")
      }
      if (length(Qvec) != length(target_vec)) {
        stop("`Qvec` and `target_vec` must have the same length.")
      }

      H_df <- private$normalize_H(H_obs, n = length(Qvec), label = "H_obs")

      is_logistic <- inherits(self$family, "family") &&
        self$family$family == "binomial" &&
        self$family$link == "logit"

      is_identity <- inherits(self$family, "family") &&
        self$family$family == "gaussian" &&
        self$family$link == "identity"

      if (is_logistic) {
        target01 <- private$clip01(self$to01(target_vec))
        offset_scale <- stats::qlogis(private$clip01(self$to01(Qvec)))
      } else {
        target01 <- target_vec
        offset_scale <- Qvec
      }

      dat <- data.frame(
        target = target01,
        offset_val = offset_scale,
        H_df,
        check.names = FALSE
      )

      H_names <- colnames(H_df)

      rhs <- if (self$use_intercept) {
        paste(H_names, collapse = " + ")
      } else {
        paste0("-1 + ", paste(H_names, collapse = " + "))
      }

      fluctuation_formula <- stats::as.formula(
        paste0("target ~ ", rhs)
      )

      fit <- suppressWarnings(stats::glm(
        formula = fluctuation_formula,
        family = self$family,
        data = dat,
        offset = dat$offset_val,
        start = rep(0, if (self$use_intercept) length(H_names) + 1L else length(H_names))
      ))

      coefs <- stats::coef(fit)

      if (self$use_intercept) {
        intercept <- unname(coefs["(Intercept)"])
        if (is.na(intercept)) intercept <- 0
        eps <- coefs[setdiff(names(coefs), "(Intercept)")]
      } else {
        intercept <- 0
        eps <- coefs
      }

      eps <- eps[H_names]
      eps[is.na(eps)] <- 0
      eps <- as.numeric(eps)
      names(eps) <- H_names

      wrap_fun <- local({
        eps_local <- eps
        intercept_local <- intercept
        H_fun_local <- H_fun_t
        Q0_fun_local <- Q0_fun_t
        is_logistic_local <- is_logistic
        clip_local <- self$clip
        to01_local <- self$to01
        from01_local <- self$from01

        function(A_vec, H_df_new) {
          Q0 <- as.numeric(Q0_fun_local(A_vec, H_df_new))
          H_new <- private$normalize_H(
            H_fun_local(A_vec, H_df_new),
            n = length(Q0),
            label = "H_fun_t(A, H)"
          )

          if (!identical(colnames(H_new), H_names)) {
            stop(
              "`H_fun_t` returned clever covariate columns that do not match ",
              "those used to fit the fluctuation."
            )
          }

          linpred_shift <- as.numeric(as.matrix(H_new) %*% eps_local) + intercept_local

          if (!is_logistic_local) {
            return(Q0 + linpred_shift)
          }

          Q001 <- pmin(pmax(to01_local(Q0), clip_local), 1 - clip_local)
          Q0star01 <- stats::plogis(stats::qlogis(Q001) + linpred_shift)
          Q0star01 <- pmin(pmax(Q0star01, clip_local), 1 - clip_local)

          from01_local(Q0star01)
        }
      })

      list(
        epsilon = eps,
        intercept = intercept,
        fit = fit,
        wrap = wrap_fun
      )
    },

    #' @description Print a compact summary.
    #' @return The object invisibly.
    print = function(...) {
      fam <- if (inherits(self$family, "family")) {
        paste0(self$family$family, "(", self$family$link, ")")
      } else {
        "<unknown>"
      }

      cat("LMTPFluctuationSubmodel\n")
      cat("  family: ", fam, "\n", sep = "")
      cat("  use_intercept: ", self$use_intercept, "\n", sep = "")
      if (fam == "binomial(logit)") {
        cat("  bounds: [", self$bounds[1], ", ", self$bounds[2], "]\n", sep = "")
        cat("  clip: ", self$clip, "\n", sep = "")
      }
      invisible(self)
    }
  ),

  private = list(
    clip01 = function(x) {
      pmin(pmax(x, self$clip), 1 - self$clip)
    },

    normalize_H = function(H, n, label = "H") {
      if (is.null(H)) {
        stop("`", label, "` cannot be NULL.")
      }

      if (is.vector(H) && !is.list(H)) {
        if (length(H) != n) {
          stop("`", label, "` must have length ", n, ".")
        }
        out <- data.frame(H = as.numeric(H))
        return(out)
      }

      if (is.matrix(H)) {
        if (nrow(H) != n) {
          stop("`", label, "` must have ", n, " rows.")
        }
        out <- as.data.frame(H, check.names = FALSE)
        if (is.null(colnames(out))) {
          colnames(out) <- paste0("H", seq_len(ncol(out)))
        }
      } else if (is.data.frame(H)) {
        if (nrow(H) != n) {
          stop("`", label, "` must have ", n, " rows.")
        }
        out <- H
        if (is.null(colnames(out))) {
          colnames(out) <- paste0("H", seq_len(ncol(out)))
        }
      } else {
        stop("`", label, "` must be a numeric vector, matrix, or data.frame.")
      }

      if (anyDuplicated(colnames(out))) {
        stop("`", label, "` must have unique column names.")
      }

      out
    }
  )
)
