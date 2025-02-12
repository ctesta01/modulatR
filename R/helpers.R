#' Check whether x is in the interval given
#'
#' `%btn%` stands for "between" and it tests whether `x` is between
#' `interval[1]` and `interval[2]`.
#'
#' The default behavior is to test if
#' $x \in [\mathtt{interval[1]}, \mathtt{interval[2]})$, but
#' to get other comparisons one can call `%btn%` without
#' using infix notation and use the `inclusivity` argument.
#'
#' @examples
#' 5 %btn% c(0, 10)
#' c(0.5, 1.5) %btn% c(0, 1)
#' `%btn%`(c(0.5, 1), c(0, 1), inclusivity = "[]")
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
#'
#' @export
`%btn[]%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "[]")
}

#' Exclusive Between
#' @export
`%btn()%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "()")
}

#' Left-Exclusive Right-Inclusive Between
#' @export
`%btn(]%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "(]")
}

#' Left-Inclusive Right-Exclusive Between
#' @export
`%btn[)%` <- function(x, interval) {
  `%btn%`(x, interval, inclusivity = "[)")
}

#' Root-Mean Squared Error
#'
#' @export
rmse <- function(x) {
  if (! is.numeric(x) || ! is.vector(x)) {
    stop("Argument x to rmse is not a numeric vector.")
  }
  return(sqrt(mean(x^2)))
}
