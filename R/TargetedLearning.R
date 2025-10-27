# Targeted Learning for an MTP

# --- One-timepoint MTP TargetParameter ----
MTPParameter <- R6::R6Class("MTPParameter", inherit = TargetParameter,
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

# --- Nuisance factory for MTP: m(A,L), m(d(A,L),L), r=g^d/g ----
MTPNuisanceFactory <- R6::R6Class("MTPNuisanceFactory",
  public = list(
    learners_Q = NULL, fml_Q = NULL,
    learners_g = NULL, fml_g = NULL,
    mtp = NULL, A_type = NULL,
    g_density_fun = NULL,
    initialize = function(learners_Q, fml_Q = NULL, fml_g = NULL,
                          mtp,
                          A_type = c("continuous","discrete"),
                          learners_g,
                          g_density_fun = NULL) {
      self$learners_Q <- learners_Q; self$fml_Q <- fml_Q
      self$learners_g <- learners_g; self$fml_g <- fml_g
      self$mtp <- mtp; self$A_type <- match.arg(A_type)
      self$g_density_fun <- g_density_fun
    },
    train = function(data, nodes) {
      # determine fml_Q and fml_g if they are NULL
      if (is.null(self$fml_Q)) {
        self$fml_Q <- paste0(nodes$Y, " ~ ", paste0(c(nodes$A, nodes$W), collapse = ' + '))
      }
      if (is.null(self$fml_g)) {
        self$fml_g <- as.formula(paste0(nodes$A, " ~ ", paste0(nodes$W, collapse = ' + ')))
      }

      # Q(A,L) via {nadir} super learner -> prediction closures
      sl_Q <- nadir::super_learner(
        data = data, formula = self$fml_Q,
        learners = self$learners_Q
      )
      QAW_fun <- function(newdata) sl_Q$predict(dplyr::mutate(newdata, Y = NA))
      QdW_fun <- function(newdata) {
        # build A* = d(A,L) rowwise using MTP
        Astar <- self$mtp$apply_policy(newdata[[nodes$A]], newdata[, setdiff(names(newdata), nodes$A), drop=FALSE])
        sl_Q$predict(dplyr::mutate(newdata, Y = NA, !!nodes$A := Astar))
      }

      if (self$A_type == "discrete") {
        # TODO: discrete case not implemented yet
        #
        # # g(A|L) with classification SL (e.g., logistic + rf); assumes {- possibly multi-level but binary typical -}
        # sl_g <- nadir::super_learner(
        #   data = data, formula = self$fml_g,
        #   learners = self$learners_g,
        #   outcome_type = 'binary'
        # )
        # g_prob <- function(nd) sl_g$predict(nd[, setdiff(names(nd), nodes$Y), drop = FALSE])
        #
        # # density ratio r(A,L) = g^d(A|L)/g(A|L) for discrete A
        # # here we assume binary A ∈ {0,1}; extendable to multi-class if learner outputs class probs
        # rAL_fun <- function(newdata) {
        #   Aobs <- newdata[[nodes$A]]
        #   Ldf  <- newdata[, setdiff(names(newdata), nodes$A), drop = FALSE]
        #   p1 <- g_prob(dplyr::mutate(newdata, !!nodes$A := 1))
        #   p0 <- 1 - p1
        #   # gd at observed A:
        #   # gd(1|L) = 1{d(1,L)=1}*p1 + 1{d(0,L)=1}*p0 ; gd(0|L) = 1{d(1,L)=0}*p1 + 1{d(0,L)=0}*p0
        #   d1 <- self$mtp$apply_policy(rep(1, nrow(newdata)), Ldf)
        #   d0 <- self$mtp$apply_policy(rep(0, nrow(newdata)), Ldf)
        #   gd1 <- as.numeric((d1 == 1)) * p1 + as.numeric((d0 == 1)) * p0
        #   gd0 <- as.numeric((d1 == 0)) * p1 + as.numeric((d0 == 0)) * p0
        #   gobs <- ifelse(Aobs == 1, p1, p0)
        #   gdobs <- ifelse(Aobs == 1, gd1, gd0)
        #   pmax(gdobs, 1e-12) / pmax(gobs, 1e-12)
        # }
        # list(
        #   QAW = QAW_fun, QdW = QdW_fun,
        #   rAL = rAL_fun,
        #   g = g_prob
        # )
      } else {
        # CONTINUOUS A
        sl_g <- nadir::super_learner(
          data = data, formula = self$fml_g,
          learners = self$learners_g
        )
        g_density <- function(nd) sl_g$predict(nd[, setdiff(names(nd), nodes$Y), drop = FALSE])
        g_density_AL <- function(A, L) {
          nd <- dplyr::bind_cols(!! nodes$A := A, L)
          g_density(nd)
        }

        rAL_fun <- function(newdata) {
          Aobs <- newdata[[nodes$A]] #newdata[nodes$A]
          Ldf  <- newdata[, setdiff(names(newdata), nodes$A), drop = FALSE]
          # g(A|L)
          g_obs <- g_density_AL(Aobs, Ldf)
          # g^d(A|L) via change-of-variables sum over inverse branches
          gd_obs <- self$mtp$gd_from_density(Aobs, Ldf, density_fun = g_density_AL)
          pmax(gd_obs, 1e-12) / pmax(g_obs, 1e-12)
        }

        list(
          QAW = QAW_fun, QdW = QdW_fun,
          rAL = rAL_fun,
          g_density = g_density_AL  # optional: expose
        )
      }
    }
  )
)

# --- TMLE fit method for MTP using your fluctuation submodels ----
fit_tmle_for_MTP <- function(data, bounds = NULL, maxit = 3, eps_tol = 1e-5) {
  nodes <- self$spec$nodes
  # Train nuisance closures
  nuis_preds <- self$nuis$train(data, nodes)

  preds <- list(
    QAW = function(newdata) nuis_preds$QAW(newdata),
    QdW = function(newdata) nuis_preds$QdW(newdata),
    rAL = function(newdata) nuis_preds$rAL(newdata)
  )

  init_preds <- preds

  A <- nodes$A; Yname <- nodes$Y
  iter <- 1; eps <- Inf
  while (iter <= maxit && eps > eps_tol) {
    cur_pred <- list(
      QAW = function(data) preds$QAW(data),
      QdW = function(data) preds$QdW(data),
      rAL = function(data) preds$rAL(data)
    )
    eif <- self$spec$eif_parts(cur_pred, data)

    # One targeting step: offset is link(QAW), covariate is H = r(A,L)
    up <- self$submodel$update_Q(
      preds = list(
        QAW = preds$QAW,
        Q1W = preds$QdW,   # not used internally but required by API; harmless
        Q0W = preds$QdW    # ditto
      ),
      eif = eif,
      data = data,
      A = A,
      Y = data[[Yname]],
      g1W = function(nd) preds$rAL(nd)  # we pass H through g1W arg slot; submodel uses it to build h
    )
    # Compose updated closures for Q(A,L); QdW remains a separate prediction of m(d(A,L),L)
    preds$QAW <- up$QAW
    eps <- up$epsilon
    iter <- iter + 1
  }

  # Final psi and IC (variance uses initial Q and r per standard practice)
  psi <- self$spec$psi(list(QdW = preds$QdW), data)

  eif <- {
    Q_init <- init_preds$QAW(data)
    r_init <- init_preds$rAL(data)
    plug   <- init_preds$QdW(data)
    (r_init * (data[[Yname]] - Q_init) + plug - psi)
  }
  se <- sqrt(stats::var(eif) / nrow(data))
  out <- list(
    psi = psi, se = se, ci95 = c(psi - 1.96*se, psi + 1.96*se),
    epsilon = up$epsilon,
    Q = list(QAW = preds$QAW, QdW = preds$QdW),
    r = preds$rAL, ic = eif
  )
  class(out) <- c("tmle_mtp","list"); out
}

