# --- One-timepoint MTP TargetParameter ----
MTPParameter <- R6::R6Class("MTPParameter", inherit = TargetedLearning::Parameter,
  public = list(
    mtp = NULL,               # instance of your MTP class (piecewise invertible)
    A_type = NULL,            # "continuous" or "discrete"
    initialize = function(mtp, A_type = c("continuous","discrete"),
                          nodes = list(W="W", A="A", Y="Y")) {
      A_type <- match.arg(A_type)
      self$mtp <- mtp; self$A_type <- A_type

      # needs: m(A,L), m(d(A,L),L), r(A,L)
      super$initialize(
        name = "LMTP1",
        nodes = nodes,
        needs = c("QAW","QdW","rAL")
      )
    },
    psi = function(preds, data) {
      # plug-in psi = E[ m{ d(A,L), L } ]
      mean(preds$QdW(data))
    },
    eif_parts = function(preds, data) {
      A <- data[[self$nodes$A]]
      Y <- data[[self$nodes$Y]]

      # clever covariate: H = r(A,L)
      H <- preds$rAL(data)
      resid <- Y - preds$QAW(data)
      plug  <- preds$QdW(data)
      list(K = matrix(H, ncol = 1), resid = resid, plug = plug)
    }
  )
)
