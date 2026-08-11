library(testthat)


# Test utilities ----

riemann_mass <- function(x,
                         y) {
  # Numerical integration helper for an equally spaced grid.
  dx <- unique(
    round(
      diff(x),
      12
    )
  )

  expect_equal(
    length(dx),
    1L
  )

  sum(y) * dx
}


make_reduce_5_mtp <- function() {
  # Canonical continuous policy:
  #
  #   d(A) = A - 5    if A >= 5
  #        = A        if A < 5.
  #
  # The image of the two branches overlaps on [0, 5), so values in that
  # interval have two valid preimages. This makes the policy useful for
  # testing the multi-branch change-of-variables calculation.

  can_reduce_5 <- function(A, H) {
    A >= 5
  }

  cannot_reduce <- function(A, H) {
    A < 5
  }

  d_reduce_5 <- function(A, H) {
    A - 5
  }

  d_identity <- function(A, H) {
    A
  }

  b_reduce_5 <- function(A_star, H) {
    A_star + 5
  }

  b_identity <- function(A_star, H) {
    A_star
  }

  db_reduce_5 <- function(A_star, H) {
    rep(
      1,
      length(A_star)
    )
  }

  db_identity <- function(A_star, H) {
    rep(
      1,
      length(A_star)
    )
  }

  MTP$new(
    treatment_type = "continuous",

    region_predicates = list(
      can_reduce_5,
      cannot_reduce
    ),

    policy_pieces = list(
      d_reduce_5,
      d_identity
    ),

    inverse_map_pieces = list(
      b_reduce_5,
      b_identity
    ),

    inverse_deriv_pieces = list(
      db_reduce_5,
      db_identity
    ),

    name = "reduce_5_if_feasible"
  )
}


# Constructor validation ----

test_that(
  "MTP constructor validates piece lengths",
  {
    always_true <- function(A, H) {
      rep(
        TRUE,
        length(A)
      )
    }

    identity_map <- function(A, H) {
      A
    }

    identity_inverse <- function(A_star, H) {
      A_star
    }

    identity_deriv <- function(A_star, H) {
      rep(
        1,
        length(A_star)
      )
    }

    expect_error(
      MTP$new(
        treatment_type = "continuous",

        region_predicates = list(
          always_true,
          always_true
        ),

        policy_pieces = list(
          identity_map
        ),

        inverse_map_pieces = list(
          identity_inverse,
          identity_inverse
        ),

        inverse_deriv_pieces = list(
          identity_deriv,
          identity_deriv
        )
      ),

      "same length"
    )
  }
)


test_that(
  "continuous MTP requires inverse maps and derivatives",
  {
    always_true <- function(A, H) {
      rep(
        TRUE,
        length(A)
      )
    }

    identity_map <- function(A, H) {
      A
    }

    expect_error(
      MTP$new(
        treatment_type = "continuous",

        region_predicates = list(
          always_true
        ),

        policy_pieces = list(
          identity_map
        )
      ),

      "must be supplied"
    )
  }
)


test_that(
  "valid continuous MTP can be constructed",
  {
    expect_no_error(
      make_reduce_5_mtp()
    )
  }
)


# Forward and inverse maps ----

test_that(
  "single-piece identity policy round-trips correctly",
  {
    always_identity <- function(A, H) {
      rep(
        TRUE,
        length(A)
      )
    }

    d_identity <- function(A, H) {
      A
    }

    b_identity <- function(A_star, H) {
      A_star
    }

    db_identity <- function(A_star, H) {
      rep(
        1,
        length(A_star)
      )
    }

    mtp <- MTP$new(
      treatment_type = "continuous",

      region_predicates = list(
        always_identity
      ),

      policy_pieces = list(
        d_identity
      ),

      inverse_map_pieces = list(
        b_identity
      ),

      inverse_deriv_pieces = list(
        db_identity
      )
    )

    set.seed(1)

    A <- runif(
      50,
      -10,
      10
    )

    H <- data.frame(
      x = rnorm(50)
    )

    A_star <- mtp$apply_policy(
      A,
      H
    )

    A_back <- mtp$inverse_map_pieces[[1]](
      A_star,
      H
    )

    expect_equal(
      A_star,
      A,
      tolerance = 1e-12
    )

    expect_equal(
      A_back,
      A,
      tolerance = 1e-12
    )

    expect_equal(
      mtp$inverse_deriv_pieces[[1]](
        A_star,
        H
      ),
      rep(
        1,
        length(A)
      )
    )
  }
)


test_that(
  "piecewise policy applies the correct forward map",
  {
    mtp <- make_reduce_5_mtp()

    A <- c(
      4,
      6,
      9.5,
      3
    )

    H <- data.frame(
      z = rep(
        1,
        length(A)
      )
    )

    A_star <- mtp$apply_policy(
      A,
      H
    )

    expect_equal(
      A_star,
      c(
        4,
        1,
        4.5,
        3
      )
    )
  }
)


test_that(
  "inverse branch functions return the expected candidate preimages",
  {
    mtp <- make_reduce_5_mtp()

    # Both branches are valid for A* in (0, 5):
    #
    #   identity branch: A = A*
    #   shifted branch:  A = A* + 5

    A_star <- c(
      2.5,
      1
    )

    H <- data.frame(
      dummy = rep(
        1,
        length(A_star)
      )
    )

    expect_equal(
      mtp$inverse_map_pieces[[1]](
        A_star,
        H
      ),
      c(
        7.5,
        6
      )
    )

    expect_equal(
      mtp$inverse_map_pieces[[2]](
        A_star,
        H
      ),
      c(
        2.5,
        1
      )
    )

    expect_equal(
      mtp$inverse_deriv_pieces[[1]](
        A_star,
        H
      ),
      rep(
        1,
        2
      )
    )

    expect_equal(
      mtp$inverse_deriv_pieces[[2]](
        A_star,
        H
      ),
      rep(
        1,
        2
      )
    )
  }
)


# Policy-region validation ----

test_that(
  "apply_policy errors when no policy region applies",
  {
    only_unit_interval <- function(A, H) {
      A >= 0 &
        A <= 1
    }

    identity_map <- function(A, H) {
      A
    }

    identity_inverse <- function(A_star, H) {
      A_star
    }

    identity_deriv <- function(A_star, H) {
      rep(
        1,
        length(A_star)
      )
    }

    mtp <- MTP$new(
      treatment_type = "continuous",

      region_predicates = list(
        only_unit_interval
      ),

      policy_pieces = list(
        identity_map
      ),

      inverse_map_pieces = list(
        identity_inverse
      ),

      inverse_deriv_pieces = list(
        identity_deriv
      )
    )

    A <- c(
      -1,
      0.5,
      2
    )

    H <- data.frame(
      z = rep(
        1,
        length(A)
      )
    )

    expect_error(
      mtp$apply_policy(
        A,
        H
      ),

      "do not belong to any policy region"
    )
  }
)


test_that(
  "apply_policy errors when policy regions overlap",
  {
    all_values <- function(A, H) {
      rep(
        TRUE,
        length(A)
      )
    }

    nonnegative <- function(A, H) {
      A >= 0
    }

    identity_map <- function(A, H) {
      A
    }

    identity_inverse <- function(A_star, H) {
      A_star
    }

    identity_deriv <- function(A_star, H) {
      rep(
        1,
        length(A_star)
      )
    }

    mtp <- MTP$new(
      treatment_type = "continuous",

      region_predicates = list(
        all_values,
        nonnegative
      ),

      policy_pieces = list(
        identity_map,
        identity_map
      ),

      inverse_map_pieces = list(
        identity_inverse,
        identity_inverse
      ),

      inverse_deriv_pieces = list(
        identity_deriv,
        identity_deriv
      )
    )

    A <- c(
      -1,
      1
    )

    H <- data.frame(
      z = rep(
        1,
        length(A)
      )
    )

    expect_error(
      mtp$apply_policy(
        A,
        H
      ),

      "more than one policy region"
    )
  }
)


# Continuous density transformation ----

test_that(
  "gd_from_density sums all valid inverse branches",
  {
    mtp <- make_reduce_5_mtp()

    g_density <- function(A, H) {
      stats::dnorm(
        A,
        mean = 4,
        sd = 2
      )
    }

    # Under the policy
    #
    #   A >= 5 -> A - 5
    #   A < 5  -> A,
    #
    # the post-policy density is:
    #
    #   a* < 0:
    #       g^d(a*) = g(a*)
    #
    #   0 <= a* < 5:
    #       g^d(a*) = g(a*) + g(a* + 5)
    #
    #   a* >= 5:
    #       g^d(a*) = g(a* + 5).

    A_star <- c(
      -1,
      2,
      7
    )

    H <- data.frame(
      dummy = rep(
        1,
        length(A_star)
      )
    )

    expected <- c(
      g_density(-1, H[1, , drop = FALSE]),

      g_density(2, H[2, , drop = FALSE]) +
        g_density(7, H[2, , drop = FALSE]),

      g_density(12, H[3, , drop = FALSE])
    )

    actual <- mtp$gd_from_density(
      A_vec = A_star,
      H_df = H,
      density_fun = g_density
    )

    expect_equal(
      actual,
      expected,
      tolerance = 1e-12
    )
  }
)


test_that(
  "gd_from_density integrates to one for a Normal base density",
  {
    mtp <- make_reduce_5_mtp()

    g_density <- function(A, H) {
      stats::dnorm(
        A,
        mean = 4,
        sd = 2
      )
    }

    grid <- seq(
      -6,
      16,
      by = 0.005
    )

    H <- data.frame(
      dummy = rep(
        1,
        length(grid)
      )
    )

    g_observed <- g_density(
      grid,
      H
    )

    g_shifted <- mtp$gd_from_density(
      A_vec = grid,
      H_df = H,
      density_fun = g_density
    )

    expect_equal(
      riemann_mass(
        grid,
        g_observed
      ),
      1,
      tolerance = 5e-4
    )

    expect_equal(
      riemann_mass(
        grid,
        g_shifted
      ),
      1,
      tolerance = 5e-4
    )
  }
)


test_that(
  "gd_from_density integrates to one for a Lognormal base density",
  {
    mtp <- make_reduce_5_mtp()

    g_density <- function(A, H) {
      stats::dlnorm(
        A,
        meanlog = 1,
        sdlog = 0.5
      )
    }

    grid <- seq(
      0,
      40,
      by = 0.002
    )

    H <- data.frame(
      dummy = rep(
        1,
        length(grid)
      )
    )

    g_observed <- g_density(
      grid,
      H
    )

    g_shifted <- mtp$gd_from_density(
      A_vec = grid,
      H_df = H,
      density_fun = g_density
    )

    expect_equal(
      riemann_mass(
        grid,
        g_observed
      ),
      1,
      tolerance = 5e-4
    )

    expect_equal(
      riemann_mass(
        grid,
        g_shifted
      ),
      1,
      tolerance = 5e-4
    )
  }
)


# Additive-shift helper ----

test_that(
  "mtp_additive_shift accepts constant treatment bounds",
  {
    # Reduce treatment by five when the resulting treatment remains >= 0.
    #
    # A = 2 cannot be shifted to -3, so it remains 2.
    # A = 7 can be shifted to 2.

    mtp <- mtp_additive_shift(
      delta = -5,
      lower = 0,
      upper = Inf
    )

    A <- c(
      2,
      7
    )

    H <- data.frame(
      dummy = rep(
        1,
        length(A)
      )
    )

    expect_equal(
      mtp$apply_policy(
        A,
        H
      ),
      c(
        2,
        2
      )
    )
  }
)


test_that(
  "mtp_additive_shift accepts history-dependent treatment bounds",
  {
    mtp <- mtp_additive_shift(
      delta = 1,

      lower = function(H) {
        H$lower
      },

      upper = function(H) {
        H$upper
      }
    )

    A <- c(
      1,
      1
    )

    H <- data.frame(
      lower = c(
        0,
        0
      ),

      upper = c(
        3,
        1.5
      )
    )

    # For observation 1, 1 + 1 = 2 is feasible.
    # For observation 2, 1 + 1 = 2 exceeds its upper bound of 1.5.

    expect_equal(
      mtp$apply_policy(
        A,
        H
      ),
      c(
        2,
        1
      )
    )
  }
)


test_that(
  "mtp_additive_shift validates bound specifications",
  {
    expect_error(
      mtp_additive_shift(
        delta = -1,
        lower = "zero"
      ),

      "numeric scalar"
    )

    mtp <- mtp_additive_shift(
      delta = -1,
      lower = 2,
      upper = 1
    )

    expect_error(
      mtp$apply_policy(
        A_vec = 3,
        H_df = data.frame(
          x = 1
        )
      ),

      "cannot exceed"
    )
  }
)


test_that(
  "mtp_additive_shift induces the expected multi-branch density",
  {
    mtp <- mtp_additive_shift(
      delta = -5,
      lower = 0
    )

    g_density <- function(A, H) {
      stats::dnorm(
        A,
        mean = 4,
        sd = 2
      )
    }

    # At A* = 2 there are two valid preimages:
    #
    #   A = 2  from the identity branch
    #   A = 7  from the shifted branch.

    H <- data.frame(
      dummy = 1
    )

    actual <- mtp$gd_from_density(
      A_vec = 2,
      H_df = H,
      density_fun = g_density
    )

    expected <-
      g_density(
        2,
        H
      ) +
      g_density(
        7,
        H
      )

    expect_equal(
      actual,
      expected,
      tolerance = 1e-12
    )
  }
)


# Discrete MTPs ----

test_that(
  "discrete MTP computes post-policy mass by summing preimages",
  {
    # Observed support and mass:
    #
    #   P(A = 0) = 0.2
    #   P(A = 1) = 0.3
    #   P(A = 2) = 0.5
    #
    # Policy:
    #
    #   0 -> 0
    #   1 -> 1
    #   2 -> 1
    #
    # Therefore:
    #
    #   P(A^d = 0) = 0.2
    #   P(A^d = 1) = 0.8
    #   P(A^d = 2) = 0.

    mtp <- mtp_discrete(
      map_fun = function(A, H) {
        ifelse(
          A == 2,
          1,
          A
        )
      },

      support = 0:2
    )

    g_mass <- function(A, H) {
      probs <- c(
        `0` = 0.2,
        `1` = 0.3,
        `2` = 0.5
      )

      unname(
        probs[
          as.character(A)
        ]
      )
    }

    A_eval <- 0:2

    H <- data.frame(
      dummy = rep(
        1,
        length(A_eval)
      )
    )

    gd <- mtp$gd_from_density(
      A_vec = A_eval,
      H_df = H,
      density_fun = g_mass
    )

    expect_equal(
      gd,
      c(
        0.2,
        0.8,
        0
      ),
      tolerance = 1e-12
    )

    expect_equal(
      sum(gd),
      1,
      tolerance = 1e-12
    )
  }
)
