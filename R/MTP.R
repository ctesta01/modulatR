MTP <- R6::R6Class("MTP",
  public = list(
    # ---- User-provided components (lists of same length)
    region_predicates      = NULL, # list of f_j(A, L) -> TRUE/FALSE, whether d(A, L) lies in I_j where d=d_j
    policy_pieces          = NULL, # list of d_j(A, L) -> A_star
    inverse_map_pieces     = NULL, # list of b_j(A_star, L) -> A     (inverse of d_j on that piece's image)
    inverse_deriv_pieces   = NULL, # list of db_j(A_star, L) -> d/dA_star b_j(A_star, L)

    # ---- Constructor
    initialize = function(region_predicates,
                          policy_pieces,
                          inverse_map_pieces,
                          inverse_deriv_pieces) {

      # Allow singletons (non-piecewise) for convenience -- turn them into a length=1 list
      if (is.function(region_predicates))    region_predicates    <- list(region_predicates)
      if (is.function(policy_pieces))        policy_pieces        <- list(policy_pieces)
      if (is.function(inverse_map_pieces))   inverse_map_pieces   <- list(inverse_map_pieces)
      if (is.function(inverse_deriv_pieces)) inverse_deriv_pieces <- list(inverse_deriv_pieces)

      # Pre-flight checks
      stops <- c(
        "region_predicates must be a non-empty list of functions" =
          !(is.list(region_predicates) &&
              length(region_predicates) > 0 &&
              all(sapply(region_predicates, is.function))),
        "policy_pieces must be a list of functions of the same length" =
          !(
            is.list(policy_pieces) &&
              length(policy_pieces) == length(region_predicates) &&
              all(sapply(policy_pieces, is.function))
          ),
        "inverse_map_pieces must be a list of functions of the same length" =
          !(
            is.list(inverse_map_pieces) &&
              length(inverse_map_pieces) == length(region_predicates) &&
              all(sapply(inverse_map_pieces, is.function))
          ),
        "inverse_deriv_pieces must be a list of functions of the same length" =
          !(
            is.list(inverse_deriv_pieces) &&
              length(inverse_deriv_pieces) == length(region_predicates) &&
              all(sapply(inverse_deriv_pieces, is.function))
          )
      )
      bad <- names(stops)[unlist(stops)]
      if (length(bad)) stop(paste(bad, collapse = "\n"))

      self$region_predicates    <- region_predicates
      self$policy_pieces        <- policy_pieces
      self$inverse_map_pieces   <- inverse_map_pieces
      self$inverse_deriv_pieces <- inverse_deriv_pieces
    },

    # ---- Helper: ensure L is a data frame with rows matching A
    .coerce_AL = function(A, L) {
      if (is.vector(L) || is.atomic(L)) L <- data.frame(L = L)
      if (!is.data.frame(L)) stop("L must be a data.frame (or a vector/atomic, which will be wrapped).")
      if (length(A) != nrow(L)) stop("A and L must have the same number of observations.")
      list(A = as.vector(A), L = L)
    },

    # ---- Which piece applies for a single (A, L) *by domain predicate on natural A*
    piece_index_for_natural = function(A, L) {
      if (length(A) != 1L) stop("Provide a single A for piece_index_for_natural().")
      truth <- vapply(self$region_predicates, function(f) isTRUE(f(A, L)), logical(1))
      idx <- which(truth)
      if (length(idx) == 0L) return(NA_integer_)
      if (length(idx) > 1L) stop(paste0("More than one region predicate applies to (A: ", A, ", L: ", L, ")"))
      idx[[1]]
    },

    # ---- Vectorized: which piece for each (A, L) by *natural* A
    piece_index_for_natural_vec = function(A, L) {
      AL <- self$.coerce_AL(A, L)
      A <- AL$A
      L <- AL$L
      vapply(seq_along(A), function(i) {
        self$piece_index_for_natural(A[i], L[i, , drop = FALSE]) },
        integer(1)
        )
    },

    # ---- Forward application d_j(A,L) -> A_star (vectorized)
    apply_policy = function(A, L) {
      AL <- self$.coerce_AL(A, L); A <- AL$A; L <- AL$L
      idx <- self$piece_index_for_natural_vec(A, L)
      vapply(seq_along(A), function(i) {
        j <- idx[i]
        if (is.na(j)) return(NA_real_)
        self$policy_pieces[[j]](A[i], L[i, , drop = FALSE])
      }, numeric(1))
    },

    # ---- Given (A_star, L), find the appropriate pieces
    # We try each piece: A_nat <- b_j(A_star,L); check predicate_j(A_nat,L) AND d_j(A_nat,L) == A_star (within tol)
    .piece_indices_from_image = function(A_star, L, tol = 1e-8) {
      which(vapply(seq_along(self$inverse_map_pieces), function(j) {
        a_nat <- self$inverse_map_pieces[[j]](A_star, L)
        ok_reg <- isTRUE(self$region_predicates[[j]](a_nat, L))
        if (!ok_reg) return(FALSE)
        a_fwd <- self$policy_pieces[[j]](a_nat, L)
        is.finite(a_nat) && is.finite(a_fwd) && abs(a_fwd - A_star) <= tol
      }, logical(1)))
    },

    # ---- Vectorized inverse: b_j(A_star, L) -> A
    apply_inverse_policy = function(A_star, L, tol = 1e-8) {
      AL <- self$.coerce_AL(A_star, L); A_star <- AL$A; L <- AL$L
      lapply(seq_along(A_star), function(i) {
        js <- self$.piece_indices_from_image(A_star[i], L[i, , drop = FALSE], tol = tol)
        if (length(js) == 0L) return(numeric(0))
        vapply(js, function(j) self$inverse_map_pieces[[j]](A_star[i], L[i, , drop = FALSE]), numeric(1))
      })
    },

    # ---- Jacobian for change-of-variables: b'_j(A_star, L)
    inverse_derivative = function(A_star, L, tol = 1e-8) {
      AL <- self$.coerce_AL(A_star, L); A_star <- AL$A; L <- AL$L
      lapply(seq_along(A_star), function(i) {
        js <- self$.piece_indices_from_image(A_star[i], L[i, , drop = FALSE], tol = tol)
        if (length(js) == 0L) return(numeric(0))
        vapply(js, function(j) self$inverse_deriv_pieces[[j]](A_star[i], L[i, , drop = FALSE]), numeric(1))
      })
    },

    # Piecewise contributions for g^d at (A*, L) given a density function g(a | L)
    # density_fun takes (a_vec, L_row_df) -> numeric vector of densities
    gd_contributions = function(A_star, L, density_fun, tol = 1e-8) {
      AL <- self$.coerce_AL(A_star, L); A_star <- AL$A; L <- AL$L
      lapply(seq_along(A_star), function(i) {
        js <- self$.piece_indices_from_image(A_star[i], L[i, , drop = FALSE], tol = tol)
        if (length(js) == 0L) return(
          data.frame(piece = integer(0), a_pre = numeric(0), jac = numeric(0),
                     dens = numeric(0), term = numeric(0))
        )
        a_pre  <- vapply(js, function(j) self$inverse_map_pieces[[j]](A_star[i], L[i, , drop = FALSE]), numeric(1))
        jac    <- vapply(js, function(j) self$inverse_deriv_pieces[[j]](A_star[i], L[i, , drop = FALSE]), numeric(1))
        dens   <- density_fun(a_pre, L[i, , drop = FALSE])
        data.frame(piece = js, a_pre = a_pre, jac = jac, dens = dens,
                   term = abs(jac) * dens, row.names = NULL)
      })
    },

    # Sum over pieces to get g^d(A* | L) (vectorized over rows)
    gd_from_density = function(A_star, L, density_fun, tol = 1e-8) {
      contribs <- self$gd_contributions(A_star, L, density_fun, tol = tol)
      vapply(contribs, function(df) if (nrow(df) == 0) 0 else sum(df$term), numeric(1))
    },

    # ---- Diagnostics: check that d_j(b_j(a*,L),L) ~= a* on a grid of a* values you pass
    validate_bijections = function(A_star_grid, L, tol = 1e-6) {
      A_back  <- self$apply_inverse_policy(A_star_grid, L)
      A_fwd   <- self$apply_policy(A_back, L)
      data.frame(A_star = A_star_grid, A_back = A_back, A_fwd = A_fwd,
                 ok = is.finite(A_back) & is.finite(A_fwd) & (abs(A_fwd - A_star_grid) <= tol))
    }
  )
)

# Helpers to coerce scalar-or-function into function(L)->scalar
#' @keywords internal
.num_to_fun <- function(x) if (is.function(x)) x else (function(L) x)

#' Additive-shift MTP Construction Helper
#'
#' d(A,L) = A + delta  if lower(L) <= A + delta <= upper(L), else A
#' Inverse branches for a* are: a* - delta (shift piece), a* (identity piece).
mtp_additive_shift <- function(delta,
                               lower = -Inf,
                               upper =  Inf) {
  lower_fun <- .num_to_fun(lower)
  upper_fun <- .num_to_fun(upper)

  # piece 1: shift feasible
  pred_shift <- function(A, L) {
    lo <- lower_fun(L); hi <- upper_fun(L)
    isTRUE(A + delta >= lo && A + delta <= hi)
  }
  d_shift   <- function(A, L) A + delta
  b_shift   <- function(A_star, L) A_star - delta   # inverse on image
  db_shift  <- function(A_star, L) 1                # d/dA* of b_shift

  # piece 2: otherwise identity
  pred_id   <- function(A, L) !isTRUE(pred_shift(A, L))
  d_id      <- function(A, L) A
  b_id      <- function(A_star, L) A_star
  db_id     <- function(A_star, L) 1

  MTP$new(
    region_predicates    = list(pred_shift, pred_id),
    policy_pieces        = list(d_shift,     d_id),
    inverse_map_pieces   = list(b_shift,     b_id),
    inverse_deriv_pieces = list(db_shift,    db_id)
  )
}

#' Multiplicative-shift MTP Construction Helper
#'
#' d(A,L) = k * A  if lower(L) <= k*A <= upper(L), else A
#' Inverse branches for a* are: a*/k (shift piece), a* (identity).
#' Note: require k != 0; for k<0 the mapping flips order but is still invertible.
mtp_multiplicative_shift <- function(k,
                                     lower = 0,
                                     upper = Inf) {
  if (k == 0) stop("k must be non-zero for invertibility.")
  lower_fun <- .num_to_fun(lower)
  upper_fun <- .num_to_fun(upper)

  pred_mult <- function(A, L) {
    lo <- lower_fun(L); hi <- upper_fun(L)
    a_star <- k * A
    isTRUE(a_star >= lo && a_star <= hi)
  }
  d_mult  <- function(A, L) k * A
  b_mult  <- function(A_star, L) A_star / k
  db_mult <- function(A_star, L) 1 / k

  pred_id <- function(A, L) !isTRUE(pred_mult(A, L))
  d_id    <- function(A, L) A
  b_id    <- function(A_star, L) A_star
  db_id   <- function(A_star, L) 1

  MTP$new(
    region_predicates    = list(pred_mult,  pred_id),
    policy_pieces        = list(d_mult,     d_id),
    inverse_map_pieces   = list(b_mult,     b_id),
    inverse_deriv_pieces = list(db_mult,    db_id)
  )
}



