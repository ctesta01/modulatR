library(here)
prefix_dir <- here("inst/simstudy")
dir.create(here("results"), showWarnings = FALSE)
dir.create(here("results/sims"), showWarnings = FALSE)
dir.create(here("results/agg"), showWarnings = FALSE)
dir.create(here("sbatch-out"), showWarnings = FALSE)

# --- configuration ---
n_vals <- c(100, 300, 500, 1000)
n_reps_total <- 1000   # total repetitions per n
delta <- -0.05         # policy shift in your example
scenario <- "gaussian_simple"  # "gaussian_simple" (analytic truth) or "tv_confounded"

# Optional: chunk into <= 10k jobs. If you need more, split by multiple arrays.
set.seed(2025)
grid <- do.call(rbind, lapply(n_vals, function(n) {
  data.frame(
    idx   = NA_integer_,    # fill below
    n     = n,
    rep   = seq_len(n_reps_total),
    seed  = sample.int(.Machine$integer.max, n_reps_total),
    delta = delta,
    scenario = scenario,
    stringsAsFactors = FALSE
  )
}))

grid$idx <- seq_len(nrow(grid))

# write the grid the array will index into
write.csv(grid, file = "results/sim_grid.csv", row.names = FALSE)
message("Wrote: results/sim_grid.csv (", nrow(grid), " rows)")
