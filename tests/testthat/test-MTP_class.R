library(testthat)

# ------------------------------------------------------------
# Tests for MTP piecewise invertible policies
# ------------------------------------------------------------

# ---------- Small utilities ----------
riemann_mass <- function(x, y) {
  # x must be equally spaced
  dx <- unique(round(diff(x), 12))
  expect_equal(length(dx), 1L)
  sum(y) * dx
}

# ---------- Canonical example policy: shift down by 5 if feasible ----------
can_reduce_5  <- function(A, L) isTRUE(A >= 5)
cannot_reduce <- function(A, L) isTRUE(A < 5)

d_reduce_5    <- function(A, L) A - 5
d_identity    <- function(A, L) A

b_reduce_5    <- function(A_star, L) A_star + 5
b_identity    <- function(A_star, L) A_star

db_reduce_5   <- function(A_star, L) 1
db_identity   <- function(A_star, L) 1

# Provide local MTP constructor if not available in the namespace already
expect_mtp_available <- function() {
  expect_true(exists("MTP"), info = "MTP R6 class must exist (with multi-branch inverse).")
}

# Optional helpers if you placed them in your package
expect_helpers_available <- function() {
  expect_true(exists("mtp_additive_shift"), info = "Helper mtp_additive_shift() should exist.")
  expect_true(exists("mtp_multiplicative_shift"), info = "Helper mtp_multiplicative_shift() should exist.")
}

# ---------- Tests ----------

test_that("Constructor checks and length matching", {
  expect_mtp_available()

  # bad lengths
  expect_error(
    MTP$new(
      region_predicates    = list(can_reduce_5, cannot_reduce),
      policy_pieces        = list(d_reduce_5),
      inverse_map_pieces   = list(b_reduce_5, b_identity),
      inverse_deriv_pieces = list(db_reduce_5, db_identity)
    ),
    "same length"
  )

  # good
  expect_no_error(
    MTP$new(
      region_predicates    = list(can_reduce_5, cannot_reduce),
      policy_pieces        = list(d_reduce_5, d_identity),
      inverse_map_pieces   = list(b_reduce_5, b_identity),
      inverse_deriv_pieces = list(db_reduce_5, db_identity)
    )
  )
})

test_that("Singleton (non-piecewise) policy round-trip works", {
  expect_mtp_available()

  always_id <- function(A, L) TRUE
  d_id <- function(A, L) A
  b_id <- function(A_star, L) A_star
  db_id <- function(A_star, L) 1

  mtp <- MTP$new(
    region_predicates    = always_id,
    policy_pieces        = d_id,
    inverse_map_pieces   = b_id,
    inverse_deriv_pieces = db_id
  )

  set.seed(1)
  A <- runif(50, -10, 10)
  L <- data.frame(x = rnorm(50))
  Astar <- mtp$apply_policy(A, L)
  Aback <- unlist(mtp$apply_inverse_policy(Astar, L))
  expect_equal(Aback, A, tolerance = 1e-10)

  # Jacobian is 1 everywhere
  Jinv <- unlist(mtp$inverse_derivative(Astar, L))
  expect_true(all(abs(Jinv - 1) < 1e-12))
})

# test_that("Piecewise policy round-trip on each piece (by-piece validation)", {
#   expect_mtp_available()
#
#   mtp <- MTP$new(
#     region_predicates    = list(can_reduce_5, cannot_reduce),
#     policy_pieces        = list(d_reduce_5,  d_identity),
#     inverse_map_pieces   = list(b_reduce_5,  b_identity),
#     inverse_deriv_pieces = list(db_reduce_5, db_identity)
#   )
#
#   # validate_bijections_by_piece() must exist on your class
#   expect_true("validate_bijections_by_piece" %in% names(mtp$.__enclos_env__$self))
#
#   A_grid_list <- list(
#     seq(5, 20, by = 0.5),
#     seq(0, 4.5, by = 0.5)
#   )
#   Ltest <- data.frame(dummy = 1)
#
#   rep <- mtp$validate_bijections_by_piece(A_grid_list, Ltest, tol = 1e-10, check_multibranch = TRUE)
#
#   # Everything on-piece should round-trip
#   expect_true(all(rep$ok_region))
#   expect_true(all(rep$ok_roundtrip))
#   expect_true(all(rep$ok_inv_round))
#   expect_true(all(rep$ok_all))
#   # and multi-branch should be >=1 (often 1 for these grids; checked later for specific a*)
# })

test_that("Multi-branch inverse returns both branches for A* in (0,5)", {
  expect_mtp_available()

  mtp <- MTP$new(
    region_predicates    = list(can_reduce_5, cannot_reduce),
    policy_pieces        = list(d_reduce_5,  d_identity),
    inverse_map_pieces   = list(b_reduce_5,  b_identity),
    inverse_deriv_pieces = list(db_reduce_5, db_identity)
  )

  # Pick A* inside (0,5)
  L <- data.frame(dummy = rep(1, 3))
  Astar <- c(2.5, 1.0, 4.9)

  inv_list <- mtp$apply_inverse_policy(Astar, L)
  # For A* in (0,5): {A*, A*+5} are valid preimages
  expect_equal(inv_list[[1]], c(7.5, 2.5), tolerance = 1e-12)
  expect_equal(inv_list[[2]], c(6.0, 1.0), tolerance = 1e-12)
  expect_equal(inv_list[[3]], c(9.9, 4.9), tolerance = 1e-12)

  # Jacobians are all ones
  jinvs <- mtp$inverse_derivative(Astar, L)
  expect_equal(unlist(jinvs), rep(1, length(unlist(inv_list))), tolerance = 1e-12)
})

test_that("gd_from_density integrates to 1 for Normal base density", {
  expect_mtp_available()

  mtp <- MTP$new(
    region_predicates    = list(can_reduce_5, cannot_reduce),
    policy_pieces        = list(d_reduce_5,  d_identity),
    inverse_map_pieces   = list(b_reduce_5,  b_identity),
    inverse_deriv_pieces = list(db_reduce_5, db_identity)
  )

  g_density <- function(a_vec, L_row) dnorm(a_vec, mean = 4, sd = 2)
  Lviz <- data.frame(dummy = 1)
  a_grid <- seq(-6, 16, by = 0.005)

  g_pre  <- g_density(a_grid, Lviz)
  g_post <- mtp$gd_from_density(a_grid, Lviz[rep(1, length(a_grid)), , drop = FALSE], g_density)

  expect_equal(riemann_mass(a_grid, g_pre),  1, tolerance = 5e-4)
  expect_equal(riemann_mass(a_grid, g_post), 1, tolerance = 5e-4)
})

test_that("gd_from_density integrates to ~1 for Lognormal base density", {
  expect_mtp_available()

  mtp <- MTP$new(
    region_predicates    = list(can_reduce_5, cannot_reduce),
    policy_pieces        = list(d_reduce_5,  d_identity),
    inverse_map_pieces   = list(b_reduce_5,  b_identity),
    inverse_deriv_pieces = list(db_reduce_5, db_identity)
  )

  # Positive support; ensure grid starts at 0
  g_density <- function(a_vec, L_row) dlnorm(a_vec, meanlog = 1, sdlog = 0.5)
  Lviz <- data.frame(dummy = 1)
  a_grid <- seq(0, 40, by = 0.002)

  g_pre  <- g_density(a_grid, Lviz)
  g_post <- mtp$gd_from_density(a_grid, Lviz[rep(1, length(a_grid)), , drop = FALSE], g_density)

  expect_equal(riemann_mass(a_grid, g_pre),  1, tolerance = 5e-4)
  expect_equal(riemann_mass(a_grid, g_post), 1, tolerance = 5e-4)
})

test_that("gd_contributions and manual summation match gd_from_density", {
  expect_mtp_available()

  mtp <- MTP$new(
    region_predicates    = list(can_reduce_5, cannot_reduce),
    policy_pieces        = list(d_reduce_5,  d_identity),
    inverse_map_pieces   = list(b_reduce_5,  b_identity),
    inverse_deriv_pieces = list(db_reduce_5, db_identity)
  )

  g_density <- function(a_vec, L_row) dnorm(a_vec, mean = 4, sd = 2)
  Lviz <- data.frame(dummy = 1)

  a_grid <- seq(-2, 8, by = 0.05)
  contribs <- mtp$gd_contributions(a_grid, Lviz[rep(1, length(a_grid)), , drop = FALSE], g_density)
  gd_sum <- vapply(contribs, function(df) if (nrow(df) == 0) 0 else sum(df$term), numeric(1))
  gd_fun <- mtp$gd_from_density(a_grid, Lviz[rep(1, length(a_grid)), , drop = FALSE], g_density)

  expect_equal(gd_sum, gd_fun, tolerance = 1e-12)
})

test_that("Vectorization: shapes and NA behavior when no piece applies", {
  expect_mtp_available()

  # Artificial policy with a hole: only A in [0,1] maps via identity; else no region
  pred <- function(A, L) isTRUE(A >= 0 && A <= 1)
  d_id <- function(A, L) A; b_id <- function(A_star, L) A_star; db_id <- function(A_star, L) 1

  mtp <- MTP$new(
    region_predicates    = list(pred),
    policy_pieces        = list(d_id),
    inverse_map_pieces   = list(b_id),
    inverse_deriv_pieces = list(db_id)
  )

  A <- c(-1, 0.5, 2)
  L <- data.frame(z = rep(1, 3))
  out <- mtp$apply_policy(A, L)
  # A=-1 or 2 should yield NA (no region predicate TRUE), middle is mapped
  expect_true(is.na(out[1]))
  expect_false(is.na(out[2]))
  expect_true(is.na(out[3]))
})

test_that("Helper constructors: additive shift with bounds", {
  expect_mtp_available(); expect_helpers_available()

  # delta = -5, bounds [0, Inf): feasible shift only when A-5 >= 0 => A >= 5
  mtp_add <- mtp_additive_shift(delta = -5, lower = 0, upper = Inf)

  A <- c(2, 7)
  L <- data.frame(dummy = rep(1, 2))
  Astar <- mtp_add$apply_policy(A, L)
  expect_equal(Astar, c(2, 2))  # identity for 2, shift for 7

  inv <- mtp_add$apply_inverse_policy(A_star = 2, L = 2)
  # For A*=2 we should get two preimages: 2 (identity) and 7 (shift piece)
  expect_true(all(sort(inv[[1]]) == c(2, 7)))

  # Density transform integrates ~1 for Normal
  g_density <- function(a_vec, L_row) dnorm(a_vec, mean = 4, sd = 2)
  grid <- seq(-6, 16, by = 0.01)
  gd <- mtp_add$gd_from_density(grid, L[rep(1, length(grid)), , drop = FALSE], g_density)
  expect_equal(riemann_mass(grid, gd), 1, tolerance = 5e-4)
})

test_that("Helper constructors: multiplicative shift with bounds and Jacobian", {
  expect_mtp_available(); expect_helpers_available()

  # 20% reduction, bounds [0, Inf)
  mtp_mult <- mtp_multiplicative_shift(k = 0.8, lower = 0, upper = Inf)

  # Check inverse derivative branch is 1/k on the shift piece
  L <- data.frame(dummy = rep(1, 2))
  # Pick A in feasible region so shift applies, and another outside feasible region
  A <- c(10, -10)
  Astar <- mtp_mult$apply_policy(A, L)        # 8
  inv_list <- mtp_mult$apply_inverse_policy(Astar, L)
  jinvs <- mtp_mult$inverse_derivative(Astar, L)

  # Multiple branches possible: identity and shift; Jacobians: 1/k (shift), 1 (identity)
  expect_true(abs(jinvs[[1]] - 1/0.8) < 1e-12)
  expect_true(abs(jinvs[[2]] - 1) < 1e-12)

  # Density transform still integrates ~1
  g_density <- function(a_vec, L_row) dlnorm(a_vec, meanlog = 1, sdlog = 0.5)
  grid <- seq(0, 60, by = 0.01)
  gd <- mtp_mult$gd_from_density(grid, L[rep(1, length(grid)), , drop = FALSE], g_density)
  expect_equal(riemann_mass(grid, gd), 1, tolerance = 5e-4)
})

test_that("gd_from_density handles empty contributions gracefully", {
  expect_mtp_available()

  # Policy with predicates that never hold (no piece); gd should be 0
  never <- function(A, L) FALSE
  d0 <- function(A, L) A
  b0 <- function(A_star, L) A_star
  db0 <- function(A_star, L) 1

  mtp <- MTP$new(
    region_predicates    = list(never),
    policy_pieces        = list(d0),
    inverse_map_pieces   = list(b0),
    inverse_deriv_pieces = list(db0)
  )

  g_density <- function(a_vec, L_row) dnorm(a_vec, 0, 1)
  x <- seq(-3, 3, by = 0.1)
  gd <- mtp$gd_from_density(x, data.frame(dummy = 1)[rep(1, length(x)), , drop = FALSE], g_density)
  expect_true(all(gd == 0))
  expect_equal(riemann_mass(x, gd), 0, tolerance = 1e-12)
})

test_that("apply_policy / apply_inverse_policy are vectorized and ordered properly", {
  expect_mtp_available()

  mtp <- MTP$new(
    region_predicates    = list(can_reduce_5, cannot_reduce),
    policy_pieces        = list(d_reduce_5,  d_identity),
    inverse_map_pieces   = list(b_reduce_5,  b_identity),
    inverse_deriv_pieces = list(db_reduce_5, db_identity)
  )

  A <- c(4, 6, 9.5, 3)
  L <- data.frame(z = rep(1, 4))
  Astar <- mtp$apply_policy(A, L)
  expect_equal(Astar, c(4, 1, 4.5, 3))

  inv <- mtp$apply_inverse_policy(Astar, L)
  # For A*=1 (second element), preimages {1,6}; for A*=4.5 (third), {5,10}
  expect_true(all(sort(inv[[2]]) == c(1, 6)))
  expect_true(all(sort(inv[[3]]) == c(4.5, 9.5)))
})
