# R/04_experiment.R

# create the experiment

experiment <- create_experiment("Overall LMTP Estimation") |>
  add_dgp(lmtp_dgp) |>
  add_method(m_modulatR)
