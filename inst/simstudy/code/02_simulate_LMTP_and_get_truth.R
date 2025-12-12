#' Forward simulator with treatment–confounder feedback and optional LMTP
#'
#' This differs from the simulate_lmtp_data_with_simcausal() function in
#' that it 1) does not use simcausaul and 2) uses a fixed t_max == 3.
#'
#' @param n          sample size
#' @param t_max      number of time points (here 3)
#' @param policy_seq either NULL (natural world) or a list of
#'                   MTP objects of length t_max (e.g. from repeat_policy_over_time)
#' @return data frame in wide format suitable for LMTP_Data_Struct
simulate_lmtp_data_under_counterfactual_policy <-
  function(n, t_max = 3, policy_seq = NULL,
           sd_for_noise = 1) {
  stopifnot(t_max == 3) # to match Y equation below

  # ---- Baseline ----
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  V  <- W2 + 1L

  # storage for covariates and treatments
  L1_mat <- matrix(NA_real_, nrow = n, ncol = t_max)
  L2_mat <- matrix(NA_real_, nrow = n, ncol = t_max)
  A_nat  <- matrix(NA_real_, nrow = n, ncol = t_max)  # natural A_t
  A_star <- matrix(NA_real_, nrow = n, ncol = t_max)  # A_t^d

  # helper to apply policy or identity
  apply_policy_t <- function(t, A_t, L1_t, L2_t, W1, W2, V) {
    if (is.null(policy_seq)) return(A_t)
    mtp_t <- policy_seq$policies[[t]]
    # Build an L data.frame passed into MTP (include anything d depends on)
    Ldf <- data.frame(
      W1 = W1,
      W2 = W2,
      V  = V,
      L1 = L1_t,
      L2 = L2_t
    )
    mtp_t$apply_policy(A_t, Ldf)
  }

  # ---- Time 1 ----
  t <- 1L
  L1_mat[, t] <- rnorm(
    n,
    mean = 0.8 * W1 + 0.5 * W2 + 0.2 * t + 0,  # no lag A at t=1
    sd   = 1
  )
  L2_mat[, t] <- rnorm(
    n,
    mean = -0.3 * W1 + 0.7 * W2 - 0.1 * t + 0,
    sd   = 1
  )
  A_nat[, t] <- rnorm(
    n,
    mean = 0.4 * W1 + 0.2 * W2 +
      0.3 * L1_mat[, t] - 0.1 * L2_mat[, t] +
      0,  # no lag A at t=1
    sd = 1
  )
  A_star[, t] <- apply_policy_t(
    t     = t,
    A_t   = A_nat[, t],
    L1_t  = L1_mat[, t],
    L2_t  = L2_mat[, t],
    W1    = W1,
    W2    = W2,
    V     = V
  )

  # ---- Times 2 and 3 ----
  for (t in 2:t_max) {
    # covariates depend on A_{t-1}^*
    L1_mat[, t] <- rnorm(
      n,
      mean = 0.8 * W1 + 0.5 * W2 + 0.2 * t +
        0.5 * A_star[, t - 1],
      sd   = 1
    )
    L2_mat[, t] <- rnorm(
      n,
      mean = -0.3 * W1 + 0.7 * W2 - 0.1 * t +
        -0.4 * A_star[, t - 1],
      sd   = 1
    )

    # natural treatment depends on current L_t and lag A^*
    A_nat[, t] <- rnorm(
      n,
      mean = 0.4 * W1 + 0.2 * W2 +
        0.3 * L1_mat[, t] - 0.1 * L2_mat[, t] +
        0.5 * A_star[, t - 1],
      sd = 1
    )

    # apply LMTP at time t
    A_star[, t] <- apply_policy_t(
      t     = t,
      A_t   = A_nat[, t],
      L1_t  = L1_mat[, t],
      L2_t  = L2_mat[, t],
      W1    = W1,
      W2    = W2,
      V     = V
    )
  }

  # ---- Outcome uses A^* (or A if policy_seq is NULL, since then A_star == A_nat) ----
  Y_mean <- (0.2 + 0.3 * W2) * rowSums(A_star) + 0.3 * W1 + 0.1 * W2
  Y <- rnorm(n, mean = Y_mean, sd = sd_for_noise)

  # ---- Wide data frame ----
  df <- data.frame(
    id = seq_len(n),
    W1 = W1,
    W2 = W2,
    V  = V,
    L11 = L1_mat[, 1], L12 = L2_mat[, 1],
    L21 = L1_mat[, 2], L22 = L2_mat[, 2],
    L31 = L1_mat[, 3], L32 = L2_mat[, 3],
    A1  = A_star[, 1],
    A2  = A_star[, 2],
    A3  = A_star[, 3],
    Y   = Y
  )

  df
}


compare_lmtp_to_natural <- function() {

  # Natural world (no policy)
  n_truth <- 2e6L # 2 million
  dat_natural <- simulate_lmtp_data_under_counterfactual_policy(
    n = n_truth, policy_seq = NULL)
  mean_Y_natural <- mean(dat_natural$Y)
  mean_Y_natural_V1 <- dat_natural |> dplyr::filter(V == 1) |> dplyr::pull(Y) |> mean()
  mean_Y_natural_V2 <- dat_natural |> dplyr::filter(V == 2) |> dplyr::pull(Y) |> mean()

  # LMTP world with shift policy
  shift_mtp   <- mtp_additive_shift(delta = -0.05)
  policy_seq  <- repeat_policy_over_time(shift_mtp, 3)

  dat_policy  <- simulate_lmtp_data_under_counterfactual_policy(
    n = n_truth, policy_seq = policy_seq)
  mean_Y_policy <- mean(dat_policy$Y)
  mean_Y_policy_V1 <- dat_policy |> dplyr::filter(V == 1) |> dplyr::pull(Y) |> mean()
  mean_Y_policy_V2 <- dat_policy |> dplyr::filter(V == 2) |> dplyr::pull(Y) |> mean()

  # Compare the results: Natural vs. Under Policy
  simulated_truth <- c(
    EY_natural = mean_Y_natural,
    EY_policy  = mean_Y_policy,
    EY_natural_V1 = mean_Y_natural_V1,
    EY_policy_V1 = mean_Y_policy_V1,
    EY_natural_V2 = mean_Y_natural_V2,
    EY_policy_V2 = mean_Y_policy_V2,
    diff       = mean_Y_policy - mean_Y_natural,
    diff_V1 = mean_Y_policy_V1 - mean_Y_natural_V1,
    diff_V2 = mean_Y_policy_V2 - mean_Y_natural_V2
  )
  simulated_truth

}


compare_sim_methods <- function(n) {

  dat1 <- simulate_lmtp_data_with_simcausal(n)
  dat2 <- simulate_lmtp_data_under_counterfactual_policy(n, policy_seq = NULL)

  c(EY1 = mean(dat1$data$Y), EY2 = mean(dat2$Y))
}
