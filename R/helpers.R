# General helpers ----------------------------------------------------------


#' Use first argument if not null, else use second argument
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


.wald_ci <- function(est,
                     se,
                     alpha = 0.05) {
  z <- stats::qnorm(
    1 - alpha / 2
  )

  c(
    lower = est - z * se,
    upper = est + z * se
  )
}


# Policy helpers -----------------------------------------------------------


#' Construct an identity treatment-policy sequence
#'
#' @param A_type Treatment type, either `"continuous"` or `"discrete"`.
#' @param discrete_support Support of the treatment when `A_type = "discrete"`.
#' @param tau Number of treatment times.
#'
#' @return An `LMTPPolicySequence`.
#' @export
identity_policy <- function(A_type = c("continuous", "discrete"),
                            discrete_support = NULL,
                            tau = 1L) {
  A_type <- match.arg(A_type)

  if (A_type == "continuous") {
    mtp_id <- MTP$new(
      treatment_type = "continuous",

      region_predicates = list(
        function(A, H) {
          rep(TRUE, length(A))
        }
      ),

      policy_pieces = list(
        function(A, H) {
          A
        }
      ),

      inverse_map_pieces = list(
        function(A_star, H) {
          A_star
        }
      ),

      inverse_deriv_pieces = list(
        function(A_star, H) {
          rep(1, length(A_star))
        }
      ),

      name = "identity"
    )

  } else {
    mtp_id <- mtp_discrete(
      map_fun = function(A, H) {
        A
      },
      support = discrete_support,
      name = "identity"
    )
  }

  repeat_policy_over_time(
    mtp_id,
    tau = tau,
    name = "identity_sequence"
  )
}


# Subgroup helpers ---------------------------------------------------------


.make_subgroup_matrix <- function(ds,
                                  subgroup_funs) {
  if (is.function(subgroup_funs)) {
    subgroup_funs <- list(
      group = subgroup_funs
    )
  }

  if (
    !is.list(subgroup_funs) ||
    !all(
      vapply(
        subgroup_funs,
        is.function,
        logical(1)
      )
    )
  ) {
    stop(
      "`subgroup_funs` must be a function or ",
      "named list of functions."
    )
  }

  out <- lapply(
    subgroup_funs,
    function(f) {
      as.numeric(
        f(ds$df)
      )
    }
  )

  out <- as.data.frame(
    out,
    check.names = FALSE
  )

  if (is.null(names(out))) {
    names(out) <- paste0(
      "group",
      seq_len(ncol(out))
    )
  }

  out
}


# Numerical helpers --------------------------------------------------------


#' Clip probabilities away from zero and one
.clip_probability <- function(x,
                              clip_probability = 1e-6) {
  pmin(
    pmax(
      x,
      clip_probability
    ),
    1 - clip_probability
  )
}


#' Truncate positive quantities away from zero
.truncate_positive <- function(x,
                               truncate_density = 1e-12) {
  pmax(
    x,
    truncate_density
  )
}


#' Truncate values to an interval
.truncate_interval <- function(x,
                               truncate = c(1e-6, Inf)) {
  pmin(
    pmax(
      x,
      truncate[1]
    ),
    truncate[2]
  )
}




# Subgroup Matrix Construction --------------------------------------------

#' Construct indicator columns for discrete conditioning variables
#'
#' Given one or more categorical baseline variables V, construct the n x J
#' indicator matrix
#'
#'   G[i, j] = I(V_i = v_j),
#'
#' where v_j is the j-th observed combination of the conditioning variables.
#'
#' For example, if
#'
#'   subgroup_cols = c("age_group", "sex"),
#'
#' columns might correspond to
#'
#'   age_group=<50, sex=Female
#'   age_group=<50, sex=Male
#'   age_group=50-65, sex=Female
#'   ...
#'
#' These groups are mutually exclusive and exhaustive over the observed data.
#'
.make_subgroup_indicators <- function(ds,
                                      subgroup_cols) {
  if (!inherits(ds, "LMTPData")) {
    stop("`ds` must inherit from `LMTPData`.")
  }

  if (!is.character(subgroup_cols) || length(subgroup_cols) < 1L) {
    stop("`subgroup_cols` must be a non-empty character vector.")
  }

  missing_cols <- setdiff(
    subgroup_cols,
    colnames(ds$df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "The following `subgroup_cols` are not present in `ds$df`: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Conditional LMTP effects are defined using baseline covariates.
  if (
    is.null(ds$W_cols) ||
    !all(subgroup_cols %in% ds$W_cols)
  ) {
    bad <- setdiff(
      subgroup_cols,
      ds$W_cols %||% character(0)
    )

    stop(
      "`subgroup_cols` must refer to baseline variables in `ds$W_cols`. ",
      "Non-baseline columns: ",
      paste(bad, collapse = ", ")
    )
  }

  V <- ds$df[
    ,
    subgroup_cols,
    drop = FALSE
  ]

  if (anyNA(V)) {
    stop(
      "`subgroup_cols` currently cannot contain missing values. ",
      "Create an explicit missing-value category first if desired."
    )
  }

  # Build a human-readable label for the observed value of V for each row.
  #
  # For one variable:
  #   age_group=<50
  #
  # For several variables:
  #   age_group=<50, sex=Female
  label_parts <- Map(
    function(x, nm) {
      paste0(
        nm,
        "=",
        as.character(x)
      )
    },
    V,
    names(V)
  )

  subgroup_label <- do.call(
    paste,
    c(
      label_parts,
      sep = ", "
    )
  )

  # Each distinct observed value/composition of V defines one subgroup.
  subgroup_names <- unique(subgroup_label)

  # Construct
  #
  #   G[i,j] = I(V_i = v_j).
  #
  G <- vapply(
    subgroup_names,
    function(g) {
      as.numeric(
        subgroup_label == g
      )
    },
    numeric(ds$n)
  )

  # vapply gives an n x J matrix, including when J = 1.
  G <- as.data.frame(
    G,
    check.names = FALSE
  )

  colnames(G) <- subgroup_names

  G
}
