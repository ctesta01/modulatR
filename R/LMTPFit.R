#' @export
LMTPFit <- R6::R6Class(
  "LMTPFit",
  public = list(
    estimate = NULL,
    var = NULL,
    se = NULL,
    ci = NULL,
    alpha = NULL,
    parameter = NULL,
    estimator = NULL,
    n = NULL,
    ic = NULL,
    eif = NULL,
    Q_init = NULL,
    Q_star = NULL,
    omega = NULL,
    eps = NULL,
    intercept = NULL,
    subgroup_names = NULL,

    initialize = function(estimate,
                          var,
                          se,
                          ci,
                          alpha = 0.05,
                          parameter = NULL,
                          estimator = "TMLE",
                          n = NULL,
                          ic = NULL,
                          eif = NULL,
                          Q_init = NULL,
                          Q_star = NULL,
                          omega = NULL,
                          eps = NULL,
                          intercept = NULL,
                          subgroup_names = NULL) {
      self$estimate <- estimate
      self$var <- var
      self$se <- se
      self$ci <- ci
      self$alpha <- alpha
      self$parameter <- parameter
      self$estimator <- estimator
      self$n <- n
      self$ic <- ic
      self$eif <- eif
      self$Q_init <- Q_init
      self$Q_star <- Q_star
      self$omega <- omega
      self$eps <- eps
      self$intercept <- intercept
      self$subgroup_names <- subgroup_names
    },

    print = function(...) {
      cat("LMTPFit\n")
      cat("  estimator: ", self$estimator, "\n", sep = "")
      cat("  parameter: ", self$parameter, "\n", sep = "")
      cat("  n: ", self$n, "\n", sep = "")

      if (length(self$estimate) == 1L) {
        cat("  estimate: ", formatC(self$estimate, digits = 4, format = "f"), "\n", sep = "")
        cat("  se: ", formatC(self$se, digits = 4, format = "f"), "\n", sep = "")
        cat(
          "  ", 100 * (1 - self$alpha), "% CI: [",
          formatC(self$ci[1], digits = 4, format = "f"), ", ",
          formatC(self$ci[2], digits = 4, format = "f"), "]\n",
          sep = ""
        )
      } else {
        nm <- self$subgroup_names %||% paste0("param", seq_along(self$estimate))
        cat("  estimates:\n")
        for (j in seq_along(self$estimate)) {
          cat(
            "    ", nm[j], ": ",
            formatC(self$estimate[j], digits = 4, format = "f"),
            " (se = ",
            formatC(self$se[j], digits = 4, format = "f"),
            ")\n",
            sep = ""
          )
        }
      }

      invisible(self)
    }
  )
)
