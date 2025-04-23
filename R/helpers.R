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

