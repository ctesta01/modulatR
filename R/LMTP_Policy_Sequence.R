# =========================
# 0) Utilities
# =========================

# small helpers
.bound01 <- function(z, eps = 1e-6) pmin(pmax(z, eps), 1 - eps)

.to01 <- function(x, a, b) (x - a) / (b - a)
.from01 <- function(z, a, b) a + (b - a) * z

# clone a list of learners per time if user supplies single spec
.rep_if_needed <- function(x, tau, repeat_bool) {

  if (repeat_bool) {
    x <- rep(list(x), tau)
  } else {
    if (length(x) != tau || ! is.list(x)) {
      stop("formulas, learner lists, and learner arguments must be length tau if not using repeat_fmls_lnrs_and_args")
    }
  }
  return(x)
}

# =========================
# 1) Policy sequencer
# =========================

LMTPPolicySequence <- R6::R6Class(
  "LMTPPolicySequence",
  public = list(
    policies = NULL,  # list of MTP objects, length tau
    initialize = function(policies) {
      stopifnot(is.list(policies), length(policies) >= 1)
      self$policies <- policies
    },
    tau = function() length(self$policies),
    apply_policy_t = function(t, A_vec, H_df) {
      # apply the t-th MTP to vector A_t with history H_t data.frame
      self$policies[[t]]$apply_policy(A_vec, H_df)
    },
    gd_from_density_t = function(t, A_vec, H_df, density_fun) {
      self$policies[[t]]$gd_from_density(A_vec, H_df, density_fun = density_fun)
    }
  )
)

# Helper: repeat the same single-time policy across tau timepoints
repeat_policy_over_time <- function(mtp, tau) LMTPPolicySequence$new(rep(list(mtp), tau))


# Generic helper: build a per-time list of *fresh* MTP instances, one for each t
# mtp_factory: function(t, A_name, H_names) -> returns a NEW MTP instance
make_per_time_policies <- function(
    exposure_cols,      # character vector: c("A1","A2",...,"AT")
    L_cols,             # list length T: each a character vector for time-varying L_t
    W_cols = NULL,      # optional vector of baseline W names
    mtp_factory         # function(t, A_name, H_names) -> MTP (new instance)
) {
  stopifnot(length(exposure_cols) == length(L_cols))
  tau <- length(exposure_cols)
  policies <- vector("list", tau)
  for (t in seq_len(tau)) {
    A_name  <- exposure_cols[[t]]
    H_names <- c(W_cols %||% character(0), unlist(L_cols[seq_len(t)], use.names = FALSE),
                 if (t > 1) exposure_cols[seq_len(t-1)] else character(0))
    # H_t = (W, L_1..L_t, A_1..A_{t-1})
    policies[[t]] <- mtp_factory(t, A_name, H_names)
  }
  LMTPPolicySequence$new(policies)
}



# =========================
# 2) Data structure (refactor)
# =========================

LMTP_Data_Struct2 <- R6::R6Class(
  "LMTP_Data_Struct2",
  public = list(
    data = NULL,
    id_col = NULL,
    n_timesteps = NULL,
    A_cols = NULL,         # character vector length T
    L_cols = NULL,         # list length T of character vectors (time-varying covariates)
    W_cols = NULL,         # baseline covariates (single character vector or NULL)
    Y_col = NULL,

    initialize = function(data, id_col, n_timesteps, A_cols, L_cols,
                          W_cols = NULL, Y_col) {
      self$data <- data
      self$id_col <- id_col
      self$n_timesteps <- as.integer(n_timesteps)
      self$A_cols <- A_cols
      self$L_cols <- L_cols
      self$W_cols <- W_cols
      self$Y_col <- Y_col
      private$validate()
    },

    tau = function() self$n_timesteps,

    # History H_t = (A_1..A_{t-1}, L_1..L_t, W) as a data.frame
    H = function(t) {
      stopifnot(t >= 1, t <= self$tau())
      A_hist <- if (t > 1) self$A_cols[seq_len(t - 1)] else character(0)
      L_hist <- unlist(self$L_cols[seq_len(t)], use.names = FALSE)
      cols <- c(A_hist, L_hist, if (!is.null(self$W_cols)) self$W_cols)
      self$data[, cols, drop = FALSE]
    },

    # Convenience: returns (A_t, H_t) as a single data.frame with A_t column named A_t
    AH = function(t) {
      df <- self$H(t)
      df[[self$A_cols[[t]]]] <- self$data[[self$A_cols[[t]]]]
      df
    },

    A = function(t) self$data[[self$A_cols[[t]]]],
    L = function(t) self$data[, self$L_cols[[t]], drop = FALSE],
    W = function() if (is.null(self$W_cols)) NULL else self$data[, self$W_cols, drop = FALSE],
    Y = function() self$data[[self$Y_col]],

    # Time-indexed view: return a list with A, H, next-H (for convenience)
    view_t = function(t) {
      list(
        A = self$A(t),
        H = self$H(t),
        AH = self$AH(t),
        H_next = if (t < self$tau()) self$H(t + 1) else NULL
      )
    }
  ),
  private = list(
    validate = function() {
      stopifnot(self$id_col %in% colnames(self$data))
      stopifnot(length(self$A_cols) == self$n_timesteps)
      stopifnot(is.list(self$L_cols), length(self$L_cols) == self$n_timesteps)
      stopifnot(all(vapply(self$L_cols, function(v) is.character(v) && length(v) >= 1, TRUE)))
      needed <- c(unlist(self$L_cols), self$A_cols, if (!is.null(self$W_cols)) self$W_cols, self$Y_col)
      if (!all(needed %in% colnames(self$data))) {
        miss <- setdiff(needed, colnames(self$data))
        stop("LMTP_Data_Struct2: missing columns in data: ", paste(miss, collapse = ", "))
      }
    }
  )
)

# =========================
# 3) Nuisance factories for LMTP
# =========================


# A tiny helper that normalizes "learner spec" into a predictor function
make_predictor <- function(data, formula, spec, outcome_type = c("continuous","binary","density"), extra_args = NULL) {
  outcome_type <- match.arg(outcome_type)
  # Case A: Super Learner spec (a list of learners)
  if (is.list(spec) && !is.function(spec)) {
    fit <- nadir::super_learner(
      data = data, formula = formula, learners = spec,
      outcome_type = outcome_type,
      extra_learner_args = extra_args
    )
    return(list(
      kind = "sl",
      fit = fit,
      predict = function(newdata) as.numeric(fit$predict(newdata))
    ))
  }

  if ("sl_lnr_name" %in% names(attributes(spec))) {
    predictfun <- do.call(
      what = spec,
      args = c(list(formula = formula, data = data), extra_args))
    return(list(predict = predictfun))
  }

  # Case C: fixed model function, e.g., stats::glm or mgcv::gam
  if (is.function(spec)) {
    fit <- spec(formula = formula, data = data)
    pfun <- function(newdata) {
      as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
    }
    return(list(kind = "fixed", fit = fit, predict = pfun))
  }
  stop("Unsupported model spec: pass a SL library (list) or a fitting function (e.g., glm).")
}

LMTPNuisanceFactory <- R6::R6Class(
  "LMTPNuisanceFactory",
  public = list(
    learners_g = NULL,  # list length T (or single spec) of {nadir} learner libraries for g_t
    fml_g = NULL,       # optional list length T of formulas for g_t
    learners_g_extra_args = NULL, # list length T or single spec of learner extra arguments
    A_type = NULL,      # "continuous" or "discrete" per time? assume global for now
    policy_seq = NULL,  # LMTPPolicySequence
    repeat_fmls_lnrs_args = NULL, # a binary indicator TRUE/FALSE to repeat arguments to length T

    initialize = function(# learners_Q,
      learners_g,
      policy_seq,
      A_type = c("continuous","discrete"),
      fml_g = NULL,
      learners_g_extra_args = NULL,
      repeat_fmls_lnrs_args = TRUE) {
      self$learners_g <- learners_g
      self$policy_seq <- policy_seq
      self$A_type <- match.arg(A_type)
      self$fml_g <- fml_g
      self$learners_g_extra_args <- learners_g_extra_args
      self$repeat_fmls_lnrs_args <- repeat_fmls_lnrs_args
    },

    # Train and return closures:
    #  - r_t(at, ht)
    train = function(ds) {
      tau <- ds$tau()
      learners_g <- .rep_if_needed(self$learners_g, tau, self$repeat_fmls_lnrs_args)
      learners_g_extra_args <- .rep_if_needed(self$learners_g_extra_args, tau, self$repeat_fmls_lnrs_args)
      fml_g <- if (is.null(self$fml_g)) vector("list", tau) else .rep_if_needed(self$fml_g, tau, self$repeat_fmls_lnrs_args)

      g_list <- vector("list", tau)
      r_list <- vector("list", tau)

      for (t in seq_len(tau)) {
        res <- local({
          t_i <- t
          AH_i <- ds$AH(t_i)
          Aname_i <- ds$A_cols[[t_i]]

          if (is.null(fml_g[[t_i]])) {
            fml_g[[t_i]] <- as.formula(
              paste0(Aname_i, " ~ ", paste0(colnames(ds$H(t_i)), collapse = " + "))
            )
          }

          sl_g_i <- make_predictor(
            data = AH_i, formula = fml_g[[t_i]],
            spec = learners_g[[t_i]],
            outcome_type = "density",
            extra_args = learners_g_extra_args[[t_i]]
          )

          # Use tidy-eval safe col creation: "{name}" := value
          g_pred_fun_i <- function(A_vec, H_df) {
            nd <- dplyr::mutate(H_df, "{Aname_i}" := A_vec)
            as.numeric(sl_g_i$predict(nd))
          }

          if (self$A_type == "continuous") {
            r_fun_i <- function(new_A, new_H) {
              g_obs <- g_pred_fun_i(new_A, new_H)
              gd_obs <- self$policy_seq$gd_from_density_t(
                t_i, new_A, new_H, density_fun = g_pred_fun_i
              )
              pmax(gd_obs, 1e-12) / pmax(g_obs, 1e-12)
            }
          } else {
            r_fun_i <- function(new_A, new_H) {
              n <- nrow(new_H)
              A0 <- rep(0, n); A1 <- rep(1, n)
              dA0 <- self$policy_seq$apply_policy_t(t_i, A0, new_H)
              dA1 <- self$policy_seq$apply_policy_t(t_i, A1, new_H)
              g1 <- g_pred_fun_i(A1, new_H); g0 <- g_pred_fun_i(A0, new_H)
              gd1 <- (dA0 == 1) * g0 + (dA1 == 1) * g1
              gd0 <- (dA0 == 0) * g0 + (dA1 == 0) * g1
              gobs  <- ifelse(new_A == 1, g1, g0)
              gdobs <- ifelse(new_A == 1, gd1, gd0)
              pmax(gdobs, 1e-12) / pmax(gobs, 1e-12)
            }
          }

          list(idx = t_i, g = g_pred_fun_i, r = r_fun_i)
        })

        g_list[[t]] <- res$g
        r_list[[t]] <- res$r
      }

      list(r = r_list, g = g_list)
    }

  )
)

# =========================
# 6) TMLE Engine method (LMTP)
# =========================


LMTPKProvider <- R6::R6Class(
  "LMTPKProvider",
  public = list(
    r_list = NULL, ds = NULL,
    initialize = function(r_list, ds) { self$r_list <- r_list; self$ds <- ds },
    K_obs = function(t) {
      n <- nrow(self$ds$data); K <- rep(1, n)
      for (k in seq_len(t)) K <- K * self$r_list[[k]](self$ds$A(k), self$ds$H(k))
      K
    },
    # evaluate K_t at arbitrary (A,H) for the t-th factor, observed for earlier ones
    K_eval = function(t, A_vec, H_df) {
      n <- nrow(H_df); K <- rep(1, n)
      for (k in seq_len(t)) {
        if (k == t) K <- K * self$r_list[[k]](A_vec, H_df)
        else        K <- K * self$r_list[[k]](self$ds$A(k), self$ds$H(k))
      }
      K
    }
  )
)


LMTPQTrainer <- R6::R6Class(
  "LMTPQTrainer",
  public = list(
    learners_Q = NULL, fml_Q = NULL, repeat_fmls = TRUE,
    initialize = function(learners_Q, fml_Q = NULL, repeat_fmls = TRUE) {
      self$learners_Q <- learners_Q; self$fml_Q <- fml_Q; self$repeat_fmls <- repeat_fmls
    },
    train_Q_t = function(ds, t, pseudo_outcome_vec) {
      Aname <- ds$A_cols[[t]]; Ht <- ds$H(t); AHt <- ds$AH(t)
      tau <- ds$tau()

      fmls <- if (is.null(self$fml_Q)) vector("list", tau) else {
        if (self$repeat_fmls) rep(list(self$fml_Q), tau) else self$fml_Q
      }
      if (is.null(fmls[[t]])) {
        rhs <- paste0(c(colnames(Ht), Aname), collapse = " + ")
        fmls[[t]] <- stats::as.formula(paste0("M_next ~ ", rhs))
      }

      data_fit <- dplyr::bind_cols(AHt, M_next = pseudo_outcome_vec)
      sl_Q <- make_predictor(data = data_fit,
                             formula = fmls[[t]],
                             spec = self$learners_Q[[t]],
                             outcome_type = 'continuous')

      # closure Q_t^0(a,h)
      function(A_vec, H_df) {
        nd <- dplyr::bind_cols(H_df, !!rlang::sym(Aname) := A_vec, M_next = NA_real_)
        as.numeric(sl_Q$predict(nd))
      }
    }
  )
)


LMTPFluctuationIdentity <- R6::R6Class(
  "LMTPFluctuationIdentity",
  public = list(
    fit_update = function(Q0_fun_t, target_vec, K_obs_t, Kprov, t) {
      # GLM: target ~ -1 + K_obs_t + offset(Q0(A_t,H_t))
      # caller must pass Q0(A_t,H_t) as the offset argument
      fit <- stats::glm(target_vec ~ -1 + K_obs_t, family = gaussian()) # we'll supply offset externally
      # Better: fit with explicit offset vector:
      # but base R's glm needs offset as vector in model.frame; we emulate by:
      eps <- tryCatch(unname(coef(fit)[["K_obs_t"]]), error = function(e) NA_real_)
      if (is.na(eps)) eps <- 0

      # Return an updater that can evaluate at arbitrary (A,H)
      updater <- function(A_vec, H_df, offset_val) {
        offset_val + eps * Kprov$K_eval(t, A_vec, H_df)
      }

      list(epsilon = eps,
           wrap = function(A_vec, H_df) {
             Q0 <- Q0_fun_t(A_vec, H_df)
             updater(A_vec, H_df, Q0)
           })
    }
  )
)


LMTPFluctuationLogit <- R6::R6Class(
  "LMTPFluctuationLogit",
  public = list(
    a = 0, b = 1, trunc = 1e-6,
    initialize = function(bounds = c(0,1), truncation_alpha = 1e-6) {
      self$a <- bounds[1]; self$b <- bounds[2]; self$trunc <- truncation_alpha
    },
    .to01 = function(x) (x - self$a)/(self$b - self$a),
    .from01 = function(z) self$a + (self$b - self$a)*z,
    .b01 = function(z) pmin(pmax(z, self$trunc), 1 - self$trunc),

    fit_update = function(Q0_fun_t, target_vec, K_obs_t, Kprov, t) {
      Y01 <- self$.to01(target_vec)
      # offset: logit( clamp( to01( Q0(A_t,H_t) ) ) ) — supply as vector:
      # We'll compute the offset in the engine to avoid duplicating data access.
      fit <- stats::glm(Y01 ~ -1 + K_obs_t, family = binomial())
      eps <- tryCatch(unname(coef(fit)[["K_obs_t"]]), error = function(e) NA_real_)
      if (is.na(eps)) eps <- 0
      list(
        epsilon = eps,
        wrap = function(A_vec, H_df) {
          Q0 <- Q0_fun_t(A_vec, H_df)
          z  <- qlogis(self$.b01(self$.to01(Q0))) + eps * Kprov$K_eval(t, A_vec, H_df)
          self$.from01(plogis(z))
        }
      )
    }
  )
)


fit_tmle_for_LMTP <- function(
    ds,
    policy_seq,
    learners_Q,                 # list length T (or repeatable) for m_t learners
    learners_g_factory,         # your existing factory that returns r_t list (or pass r_list)
    fml_Q = NULL,
    outcome_link = c("identity","logit"),
    bounds = c(0,1),
    maxit = 1,                  # outer iterations if you want (usually 1 is fine)
    eps_tol = 1e-6,
    repeat_lnrs = TRUE
) {
  outcome_link <- match.arg(outcome_link)
  tau <- ds$tau()

  # 1) Fit g and r once (can add cross-fitting)
  nuis_gr <- learners_g_factory$train(ds)               # r_list and g_list
  r_list  <- nuis_gr$r

  # 2) Helpers
  Kprov <- LMTPKProvider$new(r_list, ds)
  learners_Q <- .rep_if_needed(learners_Q, ds$tau(), repeat_lnrs)
  Qtrainer <- LMTPQTrainer$new(learners_Q = if (length(learners_Q) == 1) rep(list(learners_Q[[1]]), tau) else learners_Q,
                               fml_Q = fml_Q, repeat_fmls = is.null(fml_Q) || length(fml_Q) == 1)
  fluct <- switch(outcome_link,
                  identity = LMTPFluctuationIdentity$new(),
                  logit    = LMTPFluctuationLogit$new(bounds = bounds))

  # 3) Initialize storage for Q0 and Q* closures
  Q0_list <- vector("list", tau)
  Qstar_list <- vector("list", tau)
  Q_star <- vector("list", tau)

  # 4) Outer iterate if desired
  iter <- 1; delta <- Inf
  eps <- 1e6
  while (iter <= maxit && eps > eps_tol) {

    # (a) start downstream pseudo-outcome with Y
    tilde_next <- ds$Y()

    # (b) backward pass: t = tau,...,1
    for (t in tau:1) {
      Ht <- ds$H(t); At <- ds$A(t)

      # build *current* pseudo-outcome for regression at step t:
      # tilde_m_{t+1}(A^d_{t+1}, H_{t+1}) already includes fluctuation at t+1 (carried from previous loop iteration)
      # For t = T, this is Y (already set).
      # Train Q_t^0 on (A_t, H_t) -> tilde_next
      Q0_t <- Qtrainer$train_Q_t(ds, t, tilde_next)

      # targets for epsilon fit at observed data:
      target_vec <- if (t < tau) {
        Hnext <- ds$H(t + 1)
        Astar_next <- policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext)
        # evaluate Q_{t+1}^*(A^d_{t+1}, H_{t+1}) using the *updated* closure from previous iteration
        Q_star[[t + 1]](Astar_next, Hnext)
      } else {
        ds$Y()
      }

      # clever covariate at observed data:
      K_obs_t <- Kprov$K_obs(t)

      # fit ε_t and wrap Q_t^*:
      up <- fluct$fit_update(Q0_fun_t = Q0_t,
                             target_vec = target_vec,
                             K_obs_t = K_obs_t,
                             Kprov = Kprov, t = t)
      eps <- up$epsilon
      Q_star[[t]] <- up$wrap

      # carry upstream the *fluctuated* pseudo-outcome: tilde_m_t(a^d_t, h_t)
      Astar_t <- policy_seq$apply_policy_t(t, ds$A(t), Ht)
      tilde_next <- Q_star[[t]](Astar_t, Ht)
    }

    iter <- iter + 1
  }

  # 5) Estimand and influence curve (using initial η if you prefer; or plug-in with Q*)
  H1 <- ds$H(1)
  Astar1 <- policy_seq$apply_policy_t(1, ds$A(1), H1)
  psi_hat <- mean(Q_star[[1]](Astar1, H1))

  # EIF with initial r and unfluctuated Q0
  w <- rep(1, nrow(ds$data))
  phi <- rep(0, nrow(ds$data))
  for (t in seq_len(tau)) {
    Ht <- ds$H(t)
    mt_obs <- Q_star[[t]](ds$A(t), Ht)  # should switch to initial Q0_t ... after storing
    mnext <- if (t < tau) {
      Hnext <- ds$H(t + 1)
      Q_star[[t + 1]](policy_seq$apply_policy_t(t + 1, ds$A(t + 1), Hnext), Hnext)
    } else ds$Y()
    w <- w * r_list[[t]](ds$A(t), Ht)
    phi <- phi + w * (mnext - mt_obs)
  }
  plug <- Q_star[[1]](Astar1, H1)
  ic <- phi + plug - psi_hat
  se <- sqrt(stats::var(ic) / nrow(ds$data))

  structure(list(
    psi = psi_hat,
    se = se,
    ci95 = c(psi_hat - 1.96*se, psi_hat + 1.96*se),
    ic = ic,
    Q = Q_star,
    r = r_list
  ), class = c("tmle_lmtp","list"))
}






