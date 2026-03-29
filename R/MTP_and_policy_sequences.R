
# MTP Class Definition ----------------------------------------------------

#' Single-time modified treatment policy
#'
#' @description
#' `MTP` represents a deterministic single-time modified treatment policy.
#'
#' It supports two common cases:
#' \itemize{
#'   \item discrete treatments, where the post-policy mass function can be
#'   computed by summing over pre-images;
#'   \item continuous treatments, where the policy is assumed piecewise
#'   invertible and the post-policy density can be computed via a
#'   change-of-variables formula.
#' }
#'
#' For continuous treatments, users supply a collection of policy regions and,
#' for each region, the forward map, inverse map, and inverse derivative.
#' This is the setup used for piecewise smooth LMTPs in Díaz et al. (2023).
#'
#' @export
#' @examples
#' ## -------------------------------------------------
#' ## Continuous MTP: simple additive shift A -> A + 1
#' ## -------------------------------------------------
#'
#' n <- 100
#' H <- data.frame(W = rnorm(n))
#' A <- rnorm(n)
#'
#' mtp <- MTP$new(
#'   treatment_type = "continuous",
#'   region_predicates = list(
#'     function(A, H) rep(TRUE, length(A))
#'   ),
#'   policy_pieces = list(
#'     function(A, H) A + 1
#'   ),
#'   inverse_map_pieces = list(
#'     function(A_star, H) A_star - 1
#'   ),
#'   inverse_deriv_pieces = list(
#'     function(A_star, H) rep(1, length(A_star))
#'   ),
#'   name = "shift_plus_one"
#' )
#'
#' mtp
#'
#' # Apply policy
#' A_star <- mtp$apply_policy(A, H)
#' head(cbind(A, A_star))
#'
#' ## -------------------------------------------------
#' ## Discrete MTP: flip binary treatment
#' ## -------------------------------------------------
#'
#' A_bin <- rbinom(n, 1, 0.5)
#'
#' mtp_flip <- MTP$new(
#'   treatment_type = "discrete",
#'   region_predicates = list(
#'     function(A, H) rep(TRUE, length(A))
#'   ),
#'   policy_pieces = list(
#'     function(A, H) 1 - A
#'   ),
#'   support = c(0, 1),
#'   name = "flip_binary"
#' )
#'
#' mtp_flip
#'
#' A_star <- mtp_flip$apply_policy(A_bin, H)
#' head(cbind(A_bin, A_star))
#'
#' ## -------------------------------------------------
#' ## Discrete MTP 2: decrement non-negative integer treatment
#' ## -------------------------------------------------
#'
#' A_pois <- rpois(n, 2)
#'
#' mtp_decrement <- MTP$new(
#'   treatment_type = "discrete",
#'   region_predicates = list(
#'     function(A, H) A > 0,
#'     function(A, H) A == 0
#'   ),
#'   policy_pieces = list(
#'     function(A, H) A - 1,
#'     function(A, H) A
#'   ),
#'   name = "decrement_integer"
#' )
#'
#' mtp_decrement
#'
#' A_star <- mtp_decrement$apply_policy(A_pois, H)
#' head(cbind(A_pois, A_star))
MTP <- R6::R6Class(
  "MTP",
  public = list(
    #' @field treatment_type Either `"continuous"` or `"discrete"`.
    treatment_type = NULL,

    #' @field name Optional name for the policy.
    name = NULL,

    #' @field region_predicates List of functions `(A, H) -> logical`
    #'   defining the regions of the observed treatment space.
    region_predicates = NULL,

    #' @field policy_pieces List of functions `(A, H) -> A_star`
    #'   defining the forward map on each region.
    policy_pieces = NULL,

    #' @field inverse_map_pieces List of functions `(A_star, H) -> A`
    #'   giving the inverse map on each region.
    inverse_map_pieces = NULL,

    #' @field inverse_deriv_pieces List of functions `(A_star, H) -> numeric`
    #'   giving the derivative of the inverse map on each region.
    inverse_deriv_pieces = NULL,

    #' @field support Optional vector of support values for discrete treatment.
    support = NULL,

    #' @description Create a new single-time MTP.
    #' @param treatment_type Either `"continuous"` or `"discrete"`.
    #' @param region_predicates List of region-membership functions.
    #' @param policy_pieces List of forward-map functions.
    #' @param inverse_map_pieces List of inverse-map functions. Required for
    #'   continuous treatments.
    #' @param inverse_deriv_pieces List of inverse-derivative functions.
    #'   Required for continuous treatments.
    #' @param support Optional support for discrete treatments.
    #' @param name Optional name.
    initialize = function(treatment_type = c("continuous", "discrete"),
                          region_predicates,
                          policy_pieces,
                          inverse_map_pieces = NULL,
                          inverse_deriv_pieces = NULL,
                          support = NULL,
                          name = NULL) {
      self$treatment_type <- match.arg(treatment_type)
      self$region_predicates <- region_predicates
      self$policy_pieces <- policy_pieces
      self$inverse_map_pieces <- inverse_map_pieces
      self$inverse_deriv_pieces <- inverse_deriv_pieces
      self$support <- support
      self$name <- name %||% paste0("MTP<", self$treatment_type, ">")

      private$validate()
      invisible(self)
    },

    #' @description Apply the policy to observed treatment values.
    #' @param A_vec Vector of observed treatment values.
    #' @param H_df Data frame of histories `H_t`.
    #' @return Vector of shifted treatment values.
    apply_policy = function(A_vec, H_df) {
      private$validate_inputs(A_vec, H_df)

      region_mat <- sapply(
        seq_along(self$region_predicates),
        function(j) self$region_predicates[[j]](A_vec, H_df)
      )

      if (!is.matrix(region_mat)) {
        region_mat <- matrix(region_mat, ncol = length(self$region_predicates))
      }

      n_hits <- rowSums(region_mat)
      if (any(n_hits == 0)) {
        stop("Some observations do not belong to any policy region.")
      }
      if (any(n_hits > 1)) {
        stop("Some observations belong to more than one policy region.")
      }

      out <- vector(mode = mode(A_vec), length = length(A_vec))
      for (j in seq_along(self$policy_pieces)) {
        idx <- region_mat[, j]
        if (any(idx)) {
          out[idx] <- self$policy_pieces[[j]](A_vec[idx], H_df[idx, , drop = FALSE])
        }
      }

      out
    },
    #' @description Evaluate the density ratio r(a, h) = g^d(a | h) / g(a | h)
    #'   from an observed-data density function.
    #' @param A_vec Vector of observed treatment values.
    #' @param H_df Data frame of histories `H_t`.
    #' @param density_fun Function `(A_vec, H_df) -> g(a | h)`.
    #' @return Numeric vector of density ratios.
    ratio_from_density = function(A_vec, H_df, density_fun) {
      g_obs <- density_fun(A_vec, H_df)
      g_post <- self$gd_from_density(A_vec, H_df, density_fun = density_fun)
      g_post / g_obs
    },

    #' @description Compute the post-policy density or mass function
    #'   `g^d(a | h)` from an observed-data density or mass function.
    #' @param A_vec Vector of treatment values at which to evaluate
    #'   the post-policy distribution.
    #' @param H_df Data frame of histories `H_t`.
    #' @param density_fun Function `(A_vec, H_df) -> g(a | h)`.
    #' @return Numeric vector.
    gd_from_density = function(A_vec, H_df, density_fun) {
      private$validate_inputs(A_vec, H_df)

      if (!is.function(density_fun)) {
        stop("`density_fun` must be a function of `(A_vec, H_df)`.")
      }

      if (self$treatment_type == "continuous") {
        return(private$gd_continuous(A_vec, H_df, density_fun))
      }

      private$gd_discrete(A_vec, H_df, density_fun)
    },

    #' @description Print a compact policy summary.
    #' @return The object invisibly.
    print = function(...) {
      cat("MTP\n")
      cat("  name: ", self$name, "\n", sep = "")
      cat("  treatment_type: ", self$treatment_type, "\n", sep = "")
      cat("  n_regions: ", length(self$policy_pieces), "\n", sep = "")
      if (self$treatment_type == "discrete" && !is.null(self$support)) {
        cat("  support: ", paste(self$support, collapse = ", "), "\n", sep = "")
      }
      invisible(self)
    }
  ),

  private = list(
    validate = function() {
      if (!is.list(self$region_predicates) || length(self$region_predicates) < 1L) {
        stop("`region_predicates` must be a non-empty list.")
      }

      if (!is.list(self$policy_pieces) || length(self$policy_pieces) < 1L) {
        stop("`policy_pieces` must be a non-empty list.")
      }

      if (length(self$region_predicates) != length(self$policy_pieces)) {
        stop("`region_predicates` and `policy_pieces` must have the same length.")
      }

      ok_pred <- vapply(self$region_predicates, is.function, logical(1))
      ok_pol  <- vapply(self$policy_pieces, is.function, logical(1))

      if (!all(ok_pred)) stop("All `region_predicates` must be functions.")
      if (!all(ok_pol))  stop("All `policy_pieces` must be functions.")

      if (self$treatment_type == "continuous") {
        if (is.null(self$inverse_map_pieces) || is.null(self$inverse_deriv_pieces)) {
          stop(
            "For continuous treatments, `inverse_map_pieces` and ",
            "`inverse_deriv_pieces` must be supplied."
          )
        }

        if (!is.list(self$inverse_map_pieces) ||
            !is.list(self$inverse_deriv_pieces) ||
            length(self$inverse_map_pieces) != length(self$policy_pieces) ||
            length(self$inverse_deriv_pieces) != length(self$policy_pieces)) {
          stop(
            "For continuous treatments, inverse maps and inverse derivatives ",
            "must be lists of the same length as `policy_pieces`."
          )
        }

        ok_inv  <- vapply(self$inverse_map_pieces, is.function, logical(1))
        ok_dinv <- vapply(self$inverse_deriv_pieces, is.function, logical(1))

        if (!all(ok_inv)) {
          stop("All `inverse_map_pieces` must be functions.")
        }
        if (!all(ok_dinv)) {
          stop("All `inverse_deriv_pieces` must be functions.")
        }
      }

      if (self$treatment_type == "discrete" && is.null(self$support)) {
        warning(
          "Discrete `MTP` created without `support`. ",
          "`gd_from_density()` will require `support`."
        )
      }
    },

    validate_inputs = function(A_vec, H_df) {
      if (length(A_vec) != nrow(H_df)) {
        stop("`length(A_vec)` must equal `nrow(H_df)`.")
      }
    },

    eval_region = function(j, A_vec, H_df) {
      out <- self$region_predicates[[j]](A_vec, H_df)

      if (!is.logical(out) || length(out) != length(A_vec)) {
        stop(
          "Region predicate ", j,
          " must return a logical vector of length length(A_vec)."
        )
      }

      out
    },

    gd_continuous = function(A_vec, H_df, density_fun, tol = 1e-8) {
      out <- numeric(length(A_vec))

      for (j in seq_along(self$inverse_map_pieces)) {
        a_back <- self$inverse_map_pieces[[j]](A_vec, H_df)

        in_region <- self$region_predicates[[j]](a_back, H_df)
        if (!is.logical(in_region) || length(in_region) != length(A_vec)) {
          stop(
            "Region predicate ", j,
            " must return a logical vector of length length(A_vec)."
          )
        }

        if (!any(in_region)) {
          next
        }

        A_forward <- self$policy_pieces[[j]](
          a_back[in_region],
          H_df[in_region, , drop = FALSE]
        )

        matches_target <- abs(A_forward - A_vec[in_region]) < tol
        if (!all(matches_target)) {
          keep_idx <- which(in_region)
          in_region[keep_idx[!matches_target]] <- FALSE
        }

        if (!any(in_region)) {
          next
        }

        jac <- self$inverse_deriv_pieces[[j]](
          A_vec[in_region],
          H_df[in_region, , drop = FALSE]
        )
        g_back <- density_fun(
          a_back[in_region],
          H_df[in_region, , drop = FALSE]
        )

        if (length(g_back) != sum(in_region)) {
          stop("`density_fun` must return a vector of the same length as its input `A_vec`.")
        }
        if (length(jac) != sum(in_region)) {
          stop(
            "Inverse derivative function ", j,
            " must return a vector of the same length as its input `A_vec`."
          )
        }

        out[in_region] <- out[in_region] + g_back * abs(jac)
      }

      out
    },

    gd_discrete = function(A_vec, H_df, density_fun) {
      if (is.null(self$support)) {
        stop("For discrete treatments, `support` must be supplied.")
      }

      out <- numeric(length(A_vec))

      for (s in self$support) {
        s_vec <- rep(s, length(A_vec))
        a_star <- self$apply_policy(s_vec, H_df)
        hits <- a_star == A_vec
        if (any(hits)) {
          g_s <- density_fun(s_vec, H_df)
          out <- out + as.numeric(hits) * g_s
        }
      }

      out
    }
  )
)



# LMTPPolicySequence Class Definition -------------------------------------


#' Longitudinal modified treatment policy sequence
#'
#' @description
#' `LMTPPolicySequence` stores one single-time `MTP` for each treatment time.
#' It provides methods for applying the policy at time `t`, evaluating the
#' post-policy treatment distribution at time `t`, and checking compatibility
#' with an `LMTPData` object.
#'
#' @export
#' @examples
#' ## -------------------------------------------------
#' ## Policy sequence with two timepoints
#' ## -------------------------------------------------
#'
#' n <- 100
#'
#' df <- data.frame(
#'   W = rnorm(n),
#'   L1 = rnorm(n),
#'   A1 = rnorm(n),
#'   L2 = rnorm(n),
#'   A2 = rnorm(n),
#'   Y  = rnorm(n)
#' )
#'
#' ds <- LMTPData$new(
#'   data = df,
#'   A_cols = c("A1", "A2"),
#'   L_cols = list("L1", "L2"),
#'   W_cols = "W",
#'   Y_col = "Y"
#' )
#'
#' # Define a simple additive shift policy
#' mtp <- mtp_additive_shift(delta = 0.5)
#'
#' policy_seq <- LMTPPolicySequence$new(
#'   policies = list(mtp, mtp)
#' )
#'
#' policy_seq
#'
#' # Apply policy at time 1
#' A1_star <- policy_seq$apply_t(1, ds$A(1), ds$H(1))
#' head(cbind(A1 = ds$A(1), A1_star))
#'
#' # Apply full sequence
#' out <- policy_seq$apply_to_data(ds)
#' head(out)
LMTPPolicySequence <- R6::R6Class(
  "LMTPPolicySequence",
  public = list(
    #' @field policies List of `MTP` objects, one for each time point.
    policies = NULL,

    #' @field name Optional name for the policy sequence.
    name = NULL,

    #' @description Create a new LMTP policy sequence.
    #' @param policies List of `MTP` objects.
    #' @param name Optional name.
    initialize = function(policies, name = NULL) {
      self$policies <- policies
      self$name <- name %||% "LMTPPolicySequence"
      private$validate()
      invisible(self)
    },

    #' @description Number of time points.
    #' @return Integer `tau`.
    tau = function() {
      length(self$policies)
    },

    #' @description Return the t-th single-time policy.
    #' @param t Time index.
    #' @return An `MTP`.
    policy_t = function(t) {
      private$check_t(t)
      self$policies[[t]]
    },

    #' @description Apply the t-th policy.
    #' @param t Time index.
    #' @param A_vec Vector of observed treatment values.
    #' @param H_df Data frame of histories `H_t`.
    #' @return Vector of shifted treatment values.
    apply_t = function(t, A_vec, H_df) {
      private$check_t(t)
      self$policies[[t]]$apply_policy(A_vec, H_df)
    },

    #' @description Evaluate the post-policy density or mass at time `t`.
    #' @param t Time index.
    #' @param A_vec Vector of treatment values.
    #' @param H_df Data frame of histories `H_t`.
    #' @param density_fun Function `(A_vec, H_df) -> g_t(a | h)`.
    #' @return Numeric vector.
    gd_t = function(t, A_vec, H_df, density_fun) {
      private$check_t(t)
      self$policies[[t]]$gd_from_density(A_vec, H_df, density_fun = density_fun)
    },

    #' @description Evaluate the density ratio
    #'   `r_t(a, h) = g_t^d(a | h) / g_t(a | h)`.
    #' @param t Time index.
    #' @param A_vec Vector of treatment values.
    #' @param H_df Data frame of histories `H_t`.
    #' @param density_fun Function `(A_vec, H_df) -> g_t(a | h)`.
    #' @return Numeric vector.
    ratio_t = function(t, A_vec, H_df, density_fun) {
      private$check_t(t)
      self$policies[[t]]$ratio_from_density(A_vec, H_df, density_fun = density_fun)
    },

    #' @description Apply the full sequence to an `LMTPData` object.
    #' @param data An `LMTPData` object.
    #' @return A data.frame with id, observed treatment, and shifted treatment
    #'   at each time.
    apply_to_data = function(data) {
      if (!inherits(data, "LMTPData")) {
        stop("`data` must inherit from `LMTPData`.")
      }
      self$validate_against_data(data)

      out <- data.frame(row.names = seq_len(data$n))
      out[[data$id_col]] <- data$df[[data$id_col]]

      for (t in seq_len(self$tau())) {
        A_name <- data$A_cols[[t]]
        A_star_name <- paste0(A_name, "_star")

        out[[A_name]] <- data$A(t)
        out[[A_star_name]] <- self$apply_t(t, data$A(t), data$H(t))
      }

      out
    },

    #' @description Check that the sequence length matches an `LMTPData` object.
    #' @param data An `LMTPData` object.
    #' @return Invisibly `TRUE` if compatible.
    validate_against_data = function(data) {
      if (!inherits(data, "LMTPData")) {
        stop("`data` must inherit from `LMTPData`.")
      }

      if (self$tau() != data$tau()) {
        stop(
          "Policy sequence has tau = ", self$tau(),
          " but `LMTPData` has tau = ", data$tau(), "."
        )
      }

      invisible(TRUE)
    },

    #' @description Print a compact summary.
    #' @return The object invisibly.
    print = function(...) {
      cat("LMTPPolicySequence\n")
      cat("  name: ", self$name, "\n", sep = "")
      cat("  tau: ", self$tau(), "\n", sep = "")
      cat("  policies:\n")
      for (t in seq_len(self$tau())) {
        cat(
          "    t = ", t, ": ",
          self$policies[[t]]$name,
          " [", self$policies[[t]]$treatment_type, "]\n",
          sep = ""
        )
      }
      invisible(self)
    }
  ),

  private = list(
    validate = function() {
      if (!is.list(self$policies) || length(self$policies) < 1L) {
        stop("`policies` must be a non-empty list.")
      }

      ok <- vapply(
        self$policies,
        function(x) inherits(x, "MTP"),
        logical(1)
      )

      if (!all(ok)) {
        stop("All elements of `policies` must inherit from `MTP`.")
      }
    },

    check_t = function(t) {
      if (!is.numeric(t) || length(t) != 1L || is.na(t)) {
        stop("`t` must be a single non-missing numeric value.")
      }

      t <- as.integer(t)

      if (t < 1L || t > length(self$policies)) {
        stop("`t` must be between 1 and tau = ", length(self$policies), ".")
      }

      invisible(t)
    }
  )
)




# helpers for MTPs --------------------------------------------------------


#' Additive shift MTP for continuous treatment
#'
#' @param delta Numeric shift.
#' @param upper_fun Function of `H_df` returning the feasible upper bound.
#' @param lower_fun Optional function of `H_df` returning a lower bound.
#' @param name Optional name.
#'
#' @return An `MTP`.
#' @export
mtp_additive_shift <- function(delta,
                               upper_fun = function(H_df) rep(Inf, nrow(H_df)),
                               lower_fun = function(H_df) rep(-Inf, nrow(H_df)),
                               name = NULL) {
  pred_shift <- function(A, H_df) {
    lo <- lower_fun(H_df)
    hi <- upper_fun(H_df)
    (A + delta >= lo) & (A + delta <= hi)
  }

  pred_id <- function(A, H_df) {
    !pred_shift(A, H_df)
  }

  d_shift <- function(A, H_df) A + delta
  d_id    <- function(A, H_df) A

  b_shift  <- function(A_star, H_df) A_star - delta
  b_id     <- function(A_star, H_df) A_star

  db_shift <- function(A_star, H_df) rep(1, length(A_star))
  db_id    <- function(A_star, H_df) rep(1, length(A_star))

  MTP$new(
    treatment_type = "continuous",
    region_predicates = list(pred_shift, pred_id),
    policy_pieces = list(d_shift, d_id),
    inverse_map_pieces = list(b_shift, b_id),
    inverse_deriv_pieces = list(db_shift, db_id),
    name = name %||% paste0("additive_shift(", delta, ")")
  )
}

#' Lookup-table MTP for discrete treatment
#'
#' @param map_fun Function `(A, H_df) -> A_star`.
#' @param support Support of the observed treatment.
#' @param name Optional name.
#'
#' @return An `MTP`.
#' @export
mtp_discrete <- function(map_fun, support, name = NULL) {
  if (!is.function(map_fun)) {
    stop("`map_fun` must be a function.")
  }

  MTP$new(
    treatment_type = "discrete",
    region_predicates = list(function(A, H_df) rep(TRUE, length(A))),
    policy_pieces = list(map_fun),
    support = support,
    name = name %||% "discrete_mtp"
  )
}


# helpers for constructing LMTPPolicySequence -----------------------------


#' Repeat a single-time policy across all time points
#'
#' @param mtp An `MTP` object.
#' @param tau Number of time points.
#' @param name Optional name.
#'
#' @return An `LMTPPolicySequence`.
#' @export
repeat_policy_over_time <- function(mtp, tau, name = NULL) {
  if (!inherits(mtp, "MTP")) {
    stop("`mtp` must inherit from `MTP`.")
  }
  if (length(tau) != 1L || !is.numeric(tau) || tau < 1) {
    stop("`tau` must be a positive integer.")
  }

  LMTPPolicySequence$new(
    policies = rep(list(mtp), as.integer(tau)),
    name = name %||% paste0("repeat(", mtp$name, ")")
  )
}


