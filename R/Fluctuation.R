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
#' LMTP vector fluctuation submodel
#'
#' @description
#' `LMTPFluctuationSubmodel` implements the targeting step for vector-valued
#' LMTP sequential regressions.
#'
#' It takes untargeted predictions evaluated at the observed treatment and at
#' the policy-shifted treatment:
#'
#' - `m_obs = m_t(A_t, H_t)`
#' - `m_d   = m_t(A_t^d, H_t)`
#'
#' and returns targeted vectors:
#'
#' - `m_obs_star = m_t^*(A_t, H_t)`
#' - `m_d_star   = m_t^*(A_t^d, H_t)`
#'
#' @export
LMTPFluctuationSubmodel <- R6::R6Class(
  "LMTPFluctuationSubmodel",
  public = list(
    family = NULL,
    bounds = NULL,
    clip_probability = NULL,
    use_intercept = NULL,

    initialize = function(family = stats::gaussian(),
                          bounds = c(0, 1),
                          clip_probability = 1e-6,
                          use_intercept = FALSE) {
      is_logistic <- inherits(family, "family") &&
        family$family == "binomial" &&
        family$link == "logit"

      is_identity <- inherits(family, "family") &&
        family$family == "gaussian" &&
        family$link == "identity"

      if (!(is_identity || is_logistic)) {
        stop(
          "`LMTPFluctuationSubmodel` currently supports only ",
          "`gaussian(identity)` and `binomial(logit)`."
        )
      }

      if (!is.numeric(clip_probability) ||
          length(clip_probability) != 1L ||
          clip_probability <= 0 ||
          clip_probability >= 0.5) {
        stop("`clip_probability` must be a scalar in (0, 0.5).")
      }

      if (is_logistic) {
        if (!is.numeric(bounds) || length(bounds) != 2L) {
          stop("`bounds` must be a numeric vector of length 2.")
        }
        bounds <- sort(bounds)
      }

      self$family <- family
      self$bounds <- bounds
      self$clip_probability <- clip_probability
      self$use_intercept <- use_intercept

      invisible(self)
    },

    fit_update = function(m_obs,
                          target,
                          H_obs,
                          m_d,
                          H_d,
                          t = NULL) {
      private$validate_vector_inputs(
        m_obs = m_obs,
        target = target,
        H_obs = H_obs,
        m_d = m_d,
        H_d = H_d
      )

      H_obs_df <- private$normalize_H(H_obs, n = length(m_obs), label = "H_obs")
      H_d_df <- private$normalize_H(H_d, n = length(m_d), label = "H_d")

      if (!identical(colnames(H_obs_df), colnames(H_d_df))) {
        stop("`H_obs` and `H_d` must have the same column names.")
      }

      is_logistic <- private$is_logistic()

      if (is_logistic) {
        y <- private$to_unit_interval(target)
        offset_obs <- stats::qlogis(private$to_unit_interval(m_obs))
      } else {
        y <- target
        offset_obs <- m_obs
      }

      dat <- data.frame(
        y = y,
        offset_obs = offset_obs,
        H_obs_df,
        check.names = FALSE
      )

      H_names <- colnames(H_obs_df)

      rhs <- if (isTRUE(self$use_intercept)) {
        paste(H_names, collapse = " + ")
      } else {
        paste0("-1 + ", paste(H_names, collapse = " + "))
      }

      fluctuation_formula <- stats::as.formula(
        paste0("y ~ ", rhs)
      )

      fit <- suppressWarnings(stats::glm(
        formula = fluctuation_formula,
        family = self$family,
        data = dat,
        offset = dat$offset_obs
      ))

      coefs <- stats::coef(fit)

      if (isTRUE(self$use_intercept)) {
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

      shift_obs <- as.numeric(as.matrix(H_obs_df) %*% eps) + intercept
      shift_d <- as.numeric(as.matrix(H_d_df) %*% eps) + intercept

      if (is_logistic) {
        m_obs_star <- private$update_logistic(m_obs, shift_obs)
        m_d_star <- private$update_logistic(m_d, shift_d)
      } else {
        m_obs_star <- m_obs + shift_obs
        m_d_star <- m_d + shift_d
      }

      list(
        fit = fit,
        epsilon = eps,
        intercept = intercept,
        m_obs_star = as.numeric(m_obs_star),
        m_d_star = as.numeric(m_d_star)
      )
    },

    print = function(...) {
      fam <- paste0(self$family$family, "(", self$family$link, ")")
      cat("LMTPFluctuationSubmodel\n")
      cat("  family: ", fam, "\n", sep = "")
      cat("  use_intercept: ", self$use_intercept, "\n", sep = "")
      if (private$is_logistic()) {
        cat("  bounds: [", self$bounds[1], ", ", self$bounds[2], "]\n", sep = "")
        cat("  clip_probability: ", self$clip_probability, "\n", sep = "")
      }
      invisible(self)
    }
  ),

  private = list(
    is_logistic = function() {
      inherits(self$family, "family") &&
        self$family$family == "binomial" &&
        self$family$link == "logit"
    },

    clip01 = function(x) {
      pmin(
        pmax(x, self$clip_probability),
        1 - self$clip_probability
      )
    },

    to_unit_interval = function(x) {
      z <- (x - self$bounds[1]) / (self$bounds[2] - self$bounds[1])
      private$clip01(z)
    },

    from_unit_interval = function(z) {
      self$bounds[1] + (self$bounds[2] - self$bounds[1]) * z
    },

    update_logistic = function(m, shift) {
      m01 <- private$to_unit_interval(m)
      m_star01 <- stats::plogis(stats::qlogis(m01) + shift)
      m_star01 <- private$clip01(m_star01)
      private$from_unit_interval(m_star01)
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
      } else if (is.data.frame(H)) {
        if (nrow(H) != n) {
          stop("`", label, "` must have ", n, " rows.")
        }
        out <- H
      } else {
        stop("`", label, "` must be a numeric vector, matrix, or data.frame.")
      }

      if (is.null(colnames(out))) {
        colnames(out) <- paste0("H", seq_len(ncol(out)))
      }

      if (anyDuplicated(colnames(out))) {
        stop("`", label, "` must have unique column names.")
      }

      out
    },

    validate_vector_inputs = function(m_obs, target, H_obs, m_d, H_d) {
      if (!is.numeric(m_obs) || !is.numeric(target) || !is.numeric(m_d)) {
        stop("`m_obs`, `target`, and `m_d` must be numeric vectors.")
      }

      n <- length(m_obs)

      if (length(target) != n ||
          length(m_d) != n) {
        stop("`m_obs`, `target`, and `m_d` must have the same length.")
      }

      if (is.vector(H_obs) && !is.list(H_obs) && length(H_obs) != n) {
        stop("`H_obs` must have length ", n, ".")
      }

      if (is.vector(H_d) && !is.list(H_d) && length(H_d) != n) {
        stop("`H_d` must have length ", n, ".")
      }

      invisible(TRUE)
    }
  )
)
