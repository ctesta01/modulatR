#' Modified Treatment Policy Class Constructor
#'
#' The core idea is that the MTP class should support the storing of the
#' functions relevant to specifying a piecewise smooth invertible policy.
#'
#' In the estimation of a modified treatment policy (MTP), one needs
#' several ingredients in order to use the transformation of variables
#' implied by the policy including: the inverse policy, the segments
#' of the domain on which the policy is smooth and invertible, and the
#' derivative of the policy.
#'
#' In the initialization of an MTP, we distinguish between
#' the simplest case (no piecewise definition necessary) and the
#' more complicated case where the function is defined in a piecewise
#' manner.
#'
#' The non-piecewise case may be defined using a simplified syntax
#' that does not involve wrapping every argument inside a `list()`.
#'
#' An example of a policy where the regions must be
#'
#' @examples
#'
#' # Imagine a policy that reduces the exposure by 5 for everyone:
#'
#' mtp <- MTP$new(
#'   smooth_invertible_regions = function(A, L) { return(TRUE) },
#'   policy = function(A, L) { A - 5 },
#'   inverse_policy = function(A, L) { A + 5 },
#'   derivative_of_policy = function(A, L) { 1 },
#' )
#' mtp$policy(5)
#'
#' # An improvement grounded in reality might be to specify that the policy
#' # is a reduction of surgery duration by 5-minutes only when
#' # the original/natural surgery duration was 5+ minutes to begin with:
#'
#' mtp <- MTP$new(
#'   smooth_invertible_regions = list(
#'     region1 = function(A, L) { return( A %btn[)% c(5, Inf) ) },
#'     region2 = function(A, L) { return( A %btn[)% c(0, 5) ) }
#'     ),
#'   policy = list(
#'     policy1 = function(A, L) { A - 5 },
#'     policy2 = function(A, L) { A }
#'     ),
#'   inverse_policy = list(
#'     function(A, L) { A + 5 },
#'     function(A, L) { A }
#'   ),
#'   derivative_of_policy = list(
#'     function(A, L) { 1 },
#'     function(A, L) { 1 }
#'     )
#' )
#' mtp$policy[[2]](5)
#'
#' example_df <- data.frame(
#'   A = c(15, 4),
#'   L = c(35, 25))
#'
#' mtp$which_region(example_df$A, example_df[,'L', drop=FALSE])
#'
#' mtp$apply_policy(example_df$A, example_df[,'L', drop=FALSE])
#'
#' @exportClass MTP
MTP <- R6::R6Class("MTP",
  public = list(
    smooth_invertible_regions = NULL,
    policy = NULL,
    inverse_policy = NULL,
    derivative_of_policy = NULL,

    #' Initialize a Modified Treatment Policy (MTP)
    #'
    #' In order to specify an MTP, we require a specification of the
    #' piecewise sections the policy is invertible over (domain),
    #' the policy function (which we require to be defined piecewise
    #' over the same segments of the policy-input-space, i.e., domain),
    #' an inverse policy similarly piecewise defined, and its
    #' derivative.
    #'
    #' Here, we have opted to require that the user specify these
    #' rather than taking a computational approach to
    #' estimating them from the policy function to avoid
    #' 1) additional computing work [taking derivatives / testing for
    #' invertibility, etc. can be computationally hard] and 2) more complex code
    #' that may be challenging to debug.
    #'
    #' Defaults are not provided so that users do not forget to specify
    #' every necessary part of an MTP
    #'
    #' We don't want users to accidentally only specify part of their policy
    #' and use the defaults without realizing they've used the defaults that
    #' may not match with their policy.
    initialize = function(
      smooth_invertible_regions,
      policy,
      inverse_policy,
      derivative_of_policy
    ) {
      # Essentially we just need to store the passed objects in the
      # MTP object created.
      #
      # Additionally, we support for convenience the option to pass
      # solo functions when the policy applies equally to one region,
      # like a policy that applies "to everyone" and is smooth/invertible
      # everywhere without needing piecewise definition.   In such cases,
      # the user can just pass solo functions and the constructor will
      # automatically wrap them in an outer list() for convenience of
      # declaring simple MTPs.
      #
      if (is.function(smooth_invertible_regions)) {
        smooth_invertible_regions <- list(smooth_invertible_regions)
      }
      self$smooth_invertible_regions <- smooth_invertible_regions

      if (is.function(policy)) {
        policy <- list(policy)
      }
      self$policy <- policy

      if (is.function(inverse_policy)) {
        inverse_policy <- list(inverse_policy)
      }
      self$inverse_policy <- inverse_policy

      if (is.function(derivative_of_policy)) {
        derivative_of_policy <- list(derivative_of_policy)
      }
      self$derivative_of_policy <- derivative_of_policy
    },

  # Calculate which region the input A and L belong to
  which_region = function(A, L) {
    if (length(A) != nrow(L)) {
      stop("A and L have incompatible dimensions")
    }

    # extract the first true index
    first_true_index <- function(x) {
      true_locations <- which(x)
      if (length(true_locations) == 0) {
        return(NA)
      } else {
        return(true_locations[[1]])
      }
    }

    # for each of the test-functions that specify the smooth_invertible_regions,
    # apply it to A and L provided.
    #
    # thus for each of the functions in smooth_invertible_regions, we get a
    # TRUE/FALSE vector as long as the length of A.
    #
    # then determine the first region that applies to each A, L combination
    lapply(
      self$smooth_invertible_regions, function(f) {
        f(A, L)
      }) |>
      as.data.frame() |>
      t() |>
      apply(MARGIN = 2, first_true_index)
  },

  apply_policy = function(A, L) {

    # determine which region applies to
    which_region_applies <- self$which_region(A, L)

    # TODO: handle that we're requiring L to be a data frame of covariates
    # TODO: what happens if L is empty?
    updated_A <- sapply(
      1:length(which_region_applies),
      function(i) {
        if (is.na(which_region_applies[[i]])) {
          return(NA)
        } else {
          return(self$policy[[which_region_applies[[i]]]](A[[i]], L[i, ]))
        }
      }
      )

    return(updated_A)
  }
  )
)
