## Slated for Deprecation
##
## Not clear why this seems to produce different results from what's in the next
## file XX_simulate_LMTP_and_get_truth.R simulate_lmtp_data_under_counterfactual_policy()
##
## The problem:
##
## If one uses either to run the baseline simulation with no intervention, and
## then calculates mean(Y) at the end, one gets different results.  Downstream
## of that, one gets different results applying the LMTP methods to the data to
## estimate counterfactual means.


library(simcausal)

simulate_lmtp_data_with_simcausal <- function(n, t_max = 3L) {

  if (t_max != 3L) {
    stop("currently only t_max == 3 is supported.")
  }

  # ---- DAG specification in simcausal ----
  dag <- DAG.empty() +
    # Baseline covariates
    node("W1", t = 0, distr = "rnorm", mean = 0, sd = 1) +
    node("W2", t = 0, distr = "rbern", prob = 0.5)

  # Time-varying process with treatment–confounder feedback
  for (t in 1:t_max) {

    # L1_t depends on W1, W2, time, and A_{t-1} (for t > 1)
    dag <- add.nodes(dag, node(
      "L1", t = t, distr = "rnorm",
      mean =
        0.8 * W1[0] +
        0.5 * W2[0] +
        0.2 * t +
        ifelse(t > 1, 0.5 * A[t - 1], 0),  # A_{t-1} -> L1_t
      sd = 1
    ))

    # L2_t depends on W1, W2, time, and A_{t-1} (for t > 1)
    dag <- add.nodes(dag, node(
      "L2", t = t, distr = "rnorm",
      mean =
        -0.3 * W1[0] +
        0.7 * W2[0] -
        0.1 * t +
        ifelse(t > 1, -0.4 * A[t - 1], 0), # A_{t-1} -> L2_t
      sd = 1
    ))

    # A_t depends on W1, W2, current L1_t, L2_t, and A_{t-1}
    dag <- add.nodes(dag, node(
      "A", t = t, distr = "rnorm",
      mean =
        0.4 * W1[0] +
        0.2 * W2[0] +
        0.3 * L1[t] +
        -0.1 * L2[t] +
        ifelse(t > 1, 0.5 * A[t - 1], 0),  # persistence in A
      sd = 1
    ))
  }

  # Outcome at time t_max + 1
  # Heterogeneous effect of A by W2 (V = W2 + 1):
  #  W2 = 0 -> effect 0.2 * sum A_t
  #  W2 = 1 -> effect 0.5 * sum A_t
  dag <- add.nodes(dag, node(
    "Y", t = t_max + 1, distr = "rnorm",
    mean =
      (0.2 + 0.3 * W2[0]) * (A[1] + A[2] + A[3]) +
      0.3 * W1[0] + 0.1 * W2[0],
    sd = 1
  ))

  # Finalize and simulate
  Dset <- set.DAG(dag, verbose = FALSE)
  raw  <- sim(Dset, n = n, verbose = FALSE)

  # ---- Post-process to wide format expected by LMTP_Data_Struct ----
  df <- raw %>%
    transmute(
      id = seq_len(n),
      W1 = W1_0,
      W2 = W2_0,
      V  = W2_0 + 1L,  # V in {1,2}
      L11 = L1_1, L12 = L2_1,
      L21 = L1_2, L22 = L2_2,
      L31 = L1_3, L32 = L2_3,
      A1  = A_1,  A2  = A_2,  A3  = A_3,
      Y   = Y_4    # since t_max = 3, outcome is at t = 4
    )

  # LMTP data structure
  ds <- LMTP_Data_Struct$new(
    data       = df,
    id_col     = "id",
    n_timesteps = 3,
    A_cols     = c("A1", "A2", "A3"),
    L_cols     = list(c("L11", "L12"),
                      c("L21", "L22"),
                      c("L31", "L32")),
    W_cols     = c("W1", "V"),
    Y_col      = "Y"
  )

  ds
}

