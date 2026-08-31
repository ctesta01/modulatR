#' LMTP fluctuation submodel
#'
#' @description
#' `LMTPFluctuationSubmodel` implements one-step TMLE targeting for
#' longitudinal modified treatment policy estimators, in one of two forms.
#'
#' **Linear clever-covariate fluctuation** (`method = "linear"`, the default):
#'
#' * identity: \deqn{Q_t^*(A,H) = Q_t(A,H) + H_t(A,H)^\top \epsilon}
#' * logit:    \deqn{Q_t^*(A,H) = \mathrm{from01}[\mathrm{expit}(
#'   \mathrm{logit}(\mathrm{to01}(Q_t(A,H))) + H_t(A,H)^\top \epsilon)]}
#'
#' where \eqn{\epsilon} is fit by an (unweighted) GLM of the pseudo-outcome on
#' the clever covariate(s) H with offset Q_t.
#'
#' **Weighted-intercept fluctuation** (`method = "weighted"`): \eqn{\epsilon}
#' is fit by a GLM of the pseudo-outcome on baseline-measurable design columns
#' X (an intercept, or subgroup indicators), with case weights equal to the
#' cumulative density ratio \eqn{\omega_t} and offset Q_t:
#'
#' * identity: \deqn{Q_t^*(A,H) = Q_t(A,H) + X^\top \epsilon}
#' * logit:    analogous on the logit scale.
#'
#' Both forms solve the same efficient-score equation
#' \eqn{P_n[H_{t,j} (\tilde Y - Q_t^*)] = 0}; they differ only in the
#' finite-sample behavior. The weighted form self-normalizes by
#' \eqn{P_n[\omega]} rather than \eqn{P_n[\omega^2]} and is strongly preferred
#' when the density ratio is heavy-tailed (e.g. continuous additive shifts),
#' where the linear form carries an O(1/n) bias with a potentially enormous
#' constant under outcome-model misspecification. In weighted mode the update
#' is policy-invariant, so predictions of the ratio at shifted treatment
#' values never enter the targeting step.
#'
#' In weighted mode, callers must supply `weights` (length-n, nonnegative) and
#' may supply `design` (n x J numeric matrix of baseline-measurable columns;
#' `NULL` means a single intercept column). `use_intercept` is ignored in
#' weighted mode: the design matrix plays that role.
#'
#' This class deliberately does not implement the Wei et al. Lagrangian-dual
#' stabilization; that should live in a separate fluctuation or targeting
#' object.
#'
#' @export
LMTPFluctuationSubmodel <- R6::R6Class(
  "LMTPFluctuationSubmodel",
  public = list(
    family = NULL,
    bounds = NULL,
    clip_probability = NULL,
    use_intercept = NULL,
    method = NULL,

    initialize = function(family = stats::gaussian(),
                          bounds = c(0, 1),
                          clip_probability = 1e-6,
                          use_intercept = FALSE,
                          method = c("linear", "weighted")) {
      method <- match.arg(method)

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

      if (identical(method, "weighted") && isTRUE(use_intercept)) {
        warning(
          "`use_intercept` is ignored when `method = \"weighted\"`; ",
          "the design matrix supplies the intercept structure."
        )
      }

      self$family <- family
      self$bounds <- bounds
      self$clip_probability <- clip_probability
      self$use_intercept <- use_intercept
      self$method <- method

      invisible(self)
    },

    fit_update = function(m_obs,
                          target,
                          H_obs,
                          m_d,
                          H_d,
                          t = NULL,
                          weights = NULL,
                          design = NULL) {
      private$validate_vector_inputs(
        m_obs = m_obs,
        target = target,
        H_obs = H_obs,
        m_d = m_d,
        H_d = H_d
      )

      if (identical(self$method, "weighted")) {
        return(
          private$fit_update_weighted(
            m_obs = m_obs,
            target = target,
            m_d = m_d,
            weights = weights,
            design = design
          )
        )
      }

      private$fit_update_linear(
        m_obs = m_obs,
        target = target,
        H_obs = H_obs,
        m_d = m_d,
        H_d = H_d
      )
    },

    print = function(...) {
      fam <- paste0(self$family$family, "(", self$family$link, ")")
      cat("LMTPFluctuationSubmodel\n")
      cat("  method: ", self$method, "\n", sep = "")
      cat("  family: ", fam, "\n", sep = "")
      if (identical(self$method, "linear")) {
        cat("  use_intercept: ", self$use_intercept, "\n", sep = "")
      }
      if (private$is_logistic()) {
        cat("  bounds: [", self$bounds[1], ", ", self$bounds[2], "]\n", sep = "")
        cat("  clip_probability: ", self$clip_probability, "\n", sep = "")
      }
      invisible(self)
    }
  ),

  private = list(

    # ------------------------------------------------------------------ #
    # eps fit as a coefficient on a linear clever-covariate fluctuation  #
    # ------------------------------------------------------------------ #
    fit_update_linear = function(m_obs,
                                 target,
                                 H_obs,
                                 m_d,
                                 H_d) {
      H_obs_df <- private$normalize_H(H_obs, n = length(m_obs), label = "H_obs")
      H_d_df <- private$normalize_H(H_d, n = length(m_d), label = "H_d")

      if (!identical(colnames(H_obs_df), colnames(H_d_df))) {
        stop("`H_obs` and `H_d` must have the same column names.")
      }

      H_output_names <- colnames(H_obs_df)
      H_fit_names <- paste0("H", seq_along(H_output_names))

      colnames(H_obs_df) <- H_fit_names
      colnames(H_d_df) <- H_fit_names

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

      H_names <- H_fit_names

      rhs <- if (isTRUE(self$use_intercept)) {
        paste(H_names, collapse = " + ")
      } else {
        paste0("-1 + ", paste(H_names, collapse = " + "))
      }

      fluctuation_formula <- stats::as.formula(paste0("y ~ ", rhs))

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
      names(eps) <- H_output_names

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

    # ------------------------------------------------------------------ #
    # Weighted-intercept fluctuation.                                     #
    #                                                                     #
    # Solves, for each design column X_j,                                 #
    #                                                                     #
    #   Pn[ weights * X_j * (ytilde - Q*) ] = 0,                          #
    #                                                                     #
    # with Q* = Q + X eps (identity) or the logit analog. With            #
    # weights = omega_t and X = G (subgroup indicators, or an intercept   #
    # column), this is the same efficient-score equation the linear       #
    # fluctuation solves via the clever covariate H_j = (G_j/p_j) omega.  #
    #                                                                     #
    # The shift X eps is baseline-measurable, hence identical under the   #
    # observed and shifted treatments: no ratio predictions at shifted    #
    # inputs are used.                                                    #
    # ------------------------------------------------------------------ #
    fit_update_weighted = function(m_obs,
                                   target,
                                   m_d,
                                   weights,
                                   design) {
      n <- length(m_obs)

      if (is.null(weights)) {
        stop(
          "`method = \"weighted\"` requires `weights` ",
          "(the cumulative density-ratio vector omega_t)."
        )
      }

      w <- as.numeric(weights)

      if (length(w) != n || !all(is.finite(w))) {
        stop("`weights` must be a finite numeric vector of length ", n, ".")
      }
      if (any(w < 0)) {
        stop("`weights` must be nonnegative.")
      }
      if (sum(w) <= 0) {
        stop("`weights` must have a positive sum.")
      }

      X <- private$normalize_design(design, n = n)
      X_output_names <- colnames(X)
      X_fit_names <- paste0("X", seq_along(X_output_names))
      colnames(X) <- X_fit_names

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
        .fluct_w = w,
        as.data.frame(X, check.names = FALSE),
        check.names = FALSE
      )

      fluctuation_formula <- stats::as.formula(
        paste0("y ~ -1 + ", paste(X_fit_names, collapse = " + "))
      )

      fit <- suppressWarnings(stats::glm(
        formula = fluctuation_formula,
        family = self$family,
        data = dat,
        weights = dat$.fluct_w,
        offset = dat$offset_obs
      ))

      coefs <- stats::coef(fit)
      eps <- coefs[X_fit_names]

      if (anyNA(eps)) {
        warning(
          "Weighted fluctuation produced NA epsilon for design column(s): ",
          paste(X_output_names[is.na(eps)], collapse = ", "),
          ". These are set to 0, leaving the corresponding score(s) unsolved. ",
          "Check for empty or collinear subgroups."
        )
        eps[is.na(eps)] <- 0
      }

      eps <- as.numeric(eps)
      names(eps) <- X_output_names

      # Baseline-measurable shift: identical under observed and shifted
      # treatment.
      shift <- as.numeric(X %*% eps)

      if (is_logistic) {
        m_obs_star <- private$update_logistic(m_obs, shift)
        m_d_star <- private$update_logistic(m_d, shift)
      } else {
        m_obs_star <- m_obs + shift
        m_d_star <- m_d + shift
      }

      list(
        fit = fit,
        epsilon = eps,
        intercept = 0,
        m_obs_star = as.numeric(m_obs_star),
        m_d_star = as.numeric(m_d_star)
      )
    },

    normalize_design = function(design, n) {
      if (is.null(design)) {
        return(
          matrix(
            1,
            nrow = n,
            ncol = 1L,
            dimnames = list(NULL, "(targeting intercept)")
          )
        )
      }

      if (is.data.frame(design)) {
        design <- as.matrix(design)
      }

      if (is.vector(design) && !is.list(design)) {
        design <- matrix(as.numeric(design), ncol = 1L)
      }

      if (!is.matrix(design) || !is.numeric(design)) {
        stop("`design` must be NULL, a numeric vector, matrix, or data.frame.")
      }
      if (nrow(design) != n) {
        stop("`design` must have ", n, " rows.")
      }
      if (!all(is.finite(design))) {
        stop("`design` must contain only finite values.")
      }

      if (is.null(colnames(design))) {
        colnames(design) <- paste0("X", seq_len(ncol(design)))
      }
      if (anyDuplicated(colnames(design))) {
        stop("`design` must have unique column names.")
      }

      design
    },

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

