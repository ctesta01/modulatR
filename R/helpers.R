#' Internal utility: replicate a list or object over time
#'
#' Helper used throughout the LMTP implementation to expand a single
#' specification (e.g., a learner list or formula) to a list of length
#' `tau`, when `repeat_bool = TRUE`. If `repeat_bool = FALSE`, the input
#' must already be a list of length `tau`.
#'
#' @param x An object or list to be replicated.
#' @param tau Integer number of time points.
#' @param repeat_bool Logical; if `TRUE`, replicate `x` `tau` times; if
#'   `FALSE`, check that `x` is already a list of length `tau`.
#'
#' @return A list of length `tau`.
.rep_if_needed <- function(x, tau, repeat_bool) {

  if (repeat_bool) {
    x <- rep(list(x), tau)
  } else {
    if (length(x) != tau || !is.list(x)) {
      stop("formulas, learner lists, and learner arguments must be length tau if not using repeat_fmls_lnrs_and_args")
    }
  }
  return(x)
}



#' Check whether \code{x} is in the interval given
#'
#' `%btn%` stands for "between" and it tests whether `x` is between
#' `interval[1]` and `interval[2]`.
#'
#' The default behavior is to test if
#' `x %in% [ interval[1], interval[2] )`, but
#' to get other comparisons one can call `%btn%` without
#' using infix notation and use the `inclusivity` argument.
#'
#' @examples
#' 5 %btn% c(0, 10)
#' c(0.5, 1.5) %btn% c(0, 1)
#' `%btn%`(c(0.5, 1), c(0, 1), inclusivity = "[]")
#'
#' @param x A numeric vector of values to test against the interval.
#' @param interval An vector of two values defining an upper and lower limit.
#' @param inclusivity One of "[]", "()", "[)" or "(]" which are short-hand
#' respectively for inclusive, exclusive, left-inclusive and right-exclusive, and
#' left-exclusive and right-inclusive.
#'
#' @export
`%btn%` <- function(x, interval, inclusivity = "[)") {
  if (! is.numeric(interval) || length(interval) != 2) {
    stop("The interval passed to %btn% must be a length 2 numeric vector.")
  }

  if (! is.numeric(x)) {
    stop("The left-hand-side passed to %btn% must be a numeric vector.")
  }

  return(switch(
    inclusivity,
    "()" = x > interval[1] & x < interval[2],
    "[]" = x >= interval[1] & x <= interval[2],
    "[)" = x >= interval[1] & x < interval[2],
    "[)" = x > interval[1] & x <= interval[2]
  ))
}

#' Inclusive Between
#' @inheritParams %btn%
#' @export
`%btn[]%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "[]")
}

#' Exclusive Between
#' @inheritParams %btn%
#' @export
`%btn()%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "()")
}

#' Left-Exclusive Right-Inclusive Between
#' @inheritParams %btn%
#' @export
`%btn(]%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "(]")
}

#' Left-Inclusive Right-Exclusive Between
#' @inheritParams %btn%
#' @export
`%btn[)%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "[)")
}

#' Root-Mean Squared Error
#'
#' @param x A numeric vector to take the square, then mean, and then
#' square-root of.
#'
rmse <- function(x) {
  if (! is.numeric(x) || ! is.vector(x)) {
    stop("Argument x to rmse is not a numeric vector.")
  }
  return(sqrt(mean(x^2)))
}

mse <- function(x) {
  if (! is.numeric(x) || ! is.vector(x)) {
    stop("Argument x to rmse is not a numeric vector.")
  }
  return(mean(x^2))
}

