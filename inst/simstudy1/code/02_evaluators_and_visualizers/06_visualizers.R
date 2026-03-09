viz_coverage_by_n <- create_visualizer(
  .name = "coverage_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    # Prefer the evaluator output if present; otherwise compute on the fly.
    # dat <- eval_results
    # if (is.null(dat) || nrow(dat) == 0) {

      dat <- fit_results %>%
        group_by(.dgp_name, .method_name, sample_size) %>%
        summarize(
          coverage = mean(true_effect <= diff_ci_high & true_effect >= diff_ci_high, na.rm = TRUE),
          .groups = "drop"
        )
    # }

    plt <- ggplot(dat, aes(x = sample_size, y = coverage,
                           color = .method_name, shape = .method_name)) +
      geom_point() +
      geom_line(aes(group = .method_name)) +
      # facet_wrap(~ .dgp_name) +
      geom_hline(yintercept = 0.95, linetype = "dashed") +
      theme_bw() +
      xlab("Sample Size") +
      ylab("Coverage") +
      labs(color = "Method", shape = "Method") +
      ggtitle("Coverage of LMTP E[Y^d - Y] estimators across sample sizes and DGPs")

    print(plt)
    return(plt)
  }
)

# 1) Bias vs sample size, faceted by DGP
viz_bias_by_n <- create_visualizer(
  .name = "bias_by_n",
  .viz_fun = function(fit_results, eval_results, vary_params = NULL, ...) {

    dat <- eval_results[['lmtp_diff_eval']]

    x_breaks <- sort(unique(dat$sample_size))

    dat |>
      ggplot(aes(x = sample_size, y = bias)) +
      geom_point() +
      geom_line() + # aes(group = .method_name)) +
      # facet_wrap(~ .dgp_name, scales = 'free') +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_x_continuous(breaks = x_breaks) +
      theme_bw() +
      labs(
        x = "Sample Size",
        y = "Bias",
        color = "Method",
        shape = "Method",
        title = "Bias of LMTP Difference Estimators Across Sample Sizes"
      )
  },
  .description = "
### Bias vs sample size

Average bias \\(\\mathbb{E}[\\hat\\tau - \\tau]\\) by method and sample size, faceted by DGP.
Dashed horizontal line at 0.
",
  .doc_show = TRUE,
  .doc_options = list(width = 10, height = 6)
)


visualize_sqrtn_error_histograms <- function(fit_results, eval_results, vary_params = NULL,
                                             accuracy = .1, ...) {

  stopifnot(is.data.frame(fit_results))
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Package 'ggh4x' is required for this visualizer (facet_grid2).")
  }

  dat <- fit_results %>%
    mutate(
      sqrt_n_error = sqrt(sample_size) * (diff_psi - true_effect)
    ) %>%
    filter(!is.na(sqrt_n_error))

  # Make sample_size an ordered factor for nicer facet ordering
  ss_levels <- sort(unique(dat$sample_size))
  dat <- dat %>%
    mutate(
      sample_size = factor(
        sample_size,
        levels = ss_levels,
        labels = paste0("sample size: ", ss_levels)
      )
    )

  ggplot(dat, aes(x = sqrt_n_error, fill = .method_name)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
    geom_density(alpha = 0.25) +
    ggh4x::facet_grid2(
      .method_name ~ sample_size,
      scales = "free",
      independent = "x"
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = accuracy)) +
    theme_bw() +
    labs(
      x = expression(sqrt(n) %*% Error),
      y = "Density",
      fill = "Method/Estimator",
      title = expression(paste("Limiting distributions of ", sqrt(n) %*% Error)),
      subtitle = paste0(
        "Rows correspond to estimators; columns to increasing sample sizes."
      )
    )
}

viz_sqrtn_error_hist <- create_visualizer(
  .name = "sqrtn_error_hist",
  .viz_fun = visualize_sqrtn_error_histograms,
  .doc_options = list(width = 12, height = 8)
)

