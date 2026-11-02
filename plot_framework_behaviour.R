suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

parse_input_dir <- function(args) {
  default <- file.path(
    "output",
    "simple_framework_simulation"
  )
  if (length(args) == 0L) {
    return(default)
  }
  prefix <- "--input-dir="
  if (length(args) != 1L || !startsWith(args[[1]], prefix)) {
    stop("Use --input-dir=path")
  }
  substring(args[[1]], nchar(prefix) + 1L)
}

read_config_value <- function(config, quantity) {
  value <- config$value[config$quantity == quantity]
  if (length(value) != 1L) {
    stop("Missing configuration value for ", quantity)
  }
  as.numeric(value)
}

save_figure <- function(plot, stem, output_dir, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 200
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height
  )
}

mean_se <- function(x) {
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
}

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#444444"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

volatility_band <- function() {
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 7.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#D95F02",
    alpha = 0.09
  )
}

plot_dynamic_response <- function(adaptive, width_target, output_dir) {
  by_round <- adaptive |>
    group_by(round, high_volatility) |>
    summarise(
      predictive_width_mean = mean(predictive_width),
      predictive_width_se = sd(predictive_width) / sqrt(n()),
      total_n_mean = mean(total_n),
      total_n_se = sd(total_n) / sqrt(n()),
      posterior_width_mean = mean(posterior_width),
      .groups = "drop"
    )

  uncertainty <- ggplot(
    by_round,
    aes(x = round, y = predictive_width_mean)
  ) +
    volatility_band() +
    geom_errorbar(
      aes(
        ymin = predictive_width_mean - predictive_width_se,
        ymax = predictive_width_mean + predictive_width_se
      ),
      width = 0.16,
      color = "#176B87"
    ) +
    geom_line(color = "#176B87", linewidth = 0.85) +
    geom_point(color = "#176B87", size = 1.8) +
    scale_x_continuous(breaks = 1:12) +
    labs(
      x = "Sampling round",
      y = "Predictive interval width",
      title = "One-step predictive uncertainty",
      subtitle = "Rounds 6 and 7 have threefold innovation SD"
    )

  effort <- ggplot(by_round, aes(x = round, y = total_n_mean)) +
    volatility_band() +
    geom_errorbar(
      aes(
        ymin = total_n_mean - total_n_se,
        ymax = total_n_mean + total_n_se
      ),
      width = 0.16,
      color = "#176B87"
    ) +
    geom_line(color = "#176B87", linewidth = 0.85) +
    geom_point(color = "#176B87", size = 1.8) +
    scale_x_continuous(breaks = 1:12) +
    labs(
      x = "Sampling round",
      y = "Selected total sample size",
      title = "Adaptive sampling effort"
    )

  relationship <- ggplot(
    adaptive,
    aes(
      x = predictive_width,
      y = total_n,
      color = high_volatility
    )
  ) +
    geom_point(alpha = 0.5, size = 1.7) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.85
    ) +
    scale_color_manual(
      values = c(`FALSE` = "#176B87", `TRUE` = "#D95F02"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    labs(
      x = "Predictive interval width",
      y = "Selected total sample size",
      title = "Effort follows predictive uncertainty"
    )

  precision <- ggplot(
    adaptive,
    aes(x = round, y = posterior_width)
  ) +
    volatility_band() +
    geom_hline(
      yintercept = width_target,
      linetype = 2,
      linewidth = 0.65
    ) +
    geom_boxplot(
      aes(group = round),
      width = 0.58,
      outlier.alpha = 0.35,
      fill = "#8DB7C7",
      color = "#176B87"
    ) +
    scale_x_continuous(breaks = 1:12) +
    labs(
      x = "Sampling round",
      y = "Realized posterior width",
      title = "Posterior precision after sampling",
      subtitle = "Dashed line shows the width target"
    )

  figure <- (uncertainty / effort) | (relationship / precision)
  figure <- figure +
    plot_annotation(title = "How the adaptive policy responds over time") &
    base_theme
  save_figure(figure, "dynamic_response", output_dir, 11, 8)
  by_round
}

plot_allocation_behaviour <- function(adaptive, output_dir) {
  allocation_long <- adaptive |>
    select(replicate, round, n_A, n_B) |>
    pivot_longer(
      cols = c(n_A, n_B),
      names_to = "stratum",
      values_to = "allocation"
    ) |>
    mutate(stratum = recode(stratum, n_A = "Stratum A", n_B = "Stratum B"))

  allocation_round <- allocation_long |>
    group_by(round, stratum) |>
    summarise(
      mean_allocation = mean(allocation),
      se_allocation = sd(allocation) / sqrt(n()),
      .groups = "drop"
    )

  by_round <- ggplot(
    allocation_round,
    aes(x = round, y = mean_allocation, color = stratum)
  ) +
    volatility_band() +
    geom_errorbar(
      aes(
        ymin = mean_allocation - se_allocation,
        ymax = mean_allocation + se_allocation
      ),
      width = 0.14,
      linewidth = 0.4
    ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.7) +
    scale_color_manual(
      values = c("Stratum A" = "#176B87", "Stratum B" = "#8C3B2A")
    ) +
    labs(
      x = "Sampling round",
      y = "Mean allocated sample size",
      title = "Allocation by stratum"
    )

  allocation_diagnostic <- adaptive |>
    mutate(
      allocated_share_A = n_A / total_n,
      proxy_optimal_share_A =
        sqrt(allocation_coefficient_A) /
        (
          sqrt(allocation_coefficient_A) +
            sqrt(allocation_coefficient_B)
        )
    )

  proxy_relation <- ggplot(
    allocation_diagnostic,
    aes(
      x = proxy_optimal_share_A,
      y = allocated_share_A,
      color = total_n
    )
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "#555555"
    ) +
    geom_point(alpha = 0.65, size = 1.8) +
    scale_color_viridis_c(option = "C", end = 0.90) +
    coord_equal() +
    labs(
      x = "Continuous proxy share for stratum A",
      y = "Greedy allocated share for stratum A",
      color = "Total sample",
      title = "Variance proxy agreement",
      subtitle = "Departures reflect integer allocation and stratum floors"
    ) +
    theme(legend.position = "right")

  example <- allocation_long |>
    filter(replicate == 1)
  example_bars <- ggplot(
    example,
    aes(x = factor(round), y = allocation, fill = stratum)
  ) +
    geom_col(width = 0.72) +
    scale_fill_manual(
      values = c("Stratum A" = "#176B87", "Stratum B" = "#8C3B2A")
    ) +
    labs(
      x = "Sampling round",
      y = "Allocated sample size",
      fill = NULL,
      title = "One sequential allocation path"
    )

  figure <- by_round / (proxy_relation | example_bars)
  figure <- figure +
    plot_layout(heights = c(0.9, 1.1)) +
    plot_annotation(title = "Behaviour of the greedy allocation rule") &
    base_theme
  save_figure(figure, "greedy_allocation_behaviour", output_dir, 11, 8)
  allocation_round
}

cluster_bootstrap_calibration <- function(
  adaptive,
  bootstrap_replicates = 2000L,
  bootstrap_seed = 20260731L
) {
  calibration_data <- adaptive |>
    mutate(
      assurance_bin = cut(
        validated_assurance,
        breaks = seq(0.85, 1.01, by = 0.02),
        include.lowest = TRUE
      )
    ) |>
    filter(!is.na(assurance_bin))
  bin_levels <- levels(calibration_data$assurance_bin)
  cluster_ids <- sort(unique(calibration_data$replicate))
  cluster_summary <- calibration_data |>
    group_by(replicate, assurance_bin) |>
    summarise(
      validated_sum = sum(validated_assurance),
      successes = sum(width_success),
      decisions = n(),
      .groups = "drop"
    ) |>
    complete(
      replicate = cluster_ids,
      assurance_bin = factor(bin_levels, levels = bin_levels),
      fill = list(
        validated_sum = 0,
        successes = 0,
        decisions = 0
      )
    )

  validated_matrix <- xtabs(
    validated_sum ~ replicate + assurance_bin,
    cluster_summary,
    drop.unused.levels = FALSE
  )
  success_matrix <- xtabs(
    successes ~ replicate + assurance_bin,
    cluster_summary,
    drop.unused.levels = FALSE
  )
  decision_matrix <- xtabs(
    decisions ~ replicate + assurance_bin,
    cluster_summary,
    drop.unused.levels = FALSE
  )
  total_decisions <- colSums(decision_matrix)
  retained <- total_decisions > 0
  validated_matrix <- validated_matrix[, retained, drop = FALSE]
  success_matrix <- success_matrix[, retained, drop = FALSE]
  decision_matrix <- decision_matrix[, retained, drop = FALSE]
  total_decisions <- total_decisions[retained]

  set.seed(bootstrap_seed)
  bootstrap_success <- matrix(
    NA_real_,
    nrow = bootstrap_replicates,
    ncol = ncol(decision_matrix)
  )
  for (bootstrap_index in seq_len(bootstrap_replicates)) {
    sampled_clusters <- sample.int(
      nrow(decision_matrix),
      nrow(decision_matrix),
      replace = TRUE
    )
    denominators <- colSums(
      decision_matrix[sampled_clusters, , drop = FALSE]
    )
    numerators <- colSums(
      success_matrix[sampled_clusters, , drop = FALSE]
    )
    bootstrap_success[bootstrap_index, ] <-
      ifelse(denominators > 0, numerators / denominators, NA_real_)
  }

  data.frame(
    assurance_bin = colnames(decision_matrix),
    mean_validated_assurance =
      colSums(validated_matrix) / total_decisions,
    realized_success =
      colSums(success_matrix) / total_decisions,
    decisions = as.integer(total_decisions),
    lower = apply(
      bootstrap_success,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = apply(
      bootstrap_success,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
  )
}

plot_assurance_behaviour <- function(
  adaptive,
  assurance_curves,
  assurance_target,
  output_dir
) {
  display_rounds <- c(1, 5, 6, 7, 8, 12)
  curves_long <- assurance_curves |>
    filter(replicate == 1, round %in% display_rounds) |>
    select(
      round,
      total,
      assurance,
      sustained_assurance,
      sustained_lower_bound
    ) |>
    pivot_longer(
      cols = c(
        assurance,
        sustained_assurance,
        sustained_lower_bound
      ),
      names_to = "measure",
      values_to = "probability"
    ) |>
    mutate(
      measure = recode(
        measure,
        assurance = "Raw assurance",
        sustained_assurance = "Sustained assurance",
        sustained_lower_bound =
          "Sustained simultaneous lower bound"
      )
    )
  selected <- assurance_curves |>
    filter(replicate == 1, round %in% display_rounds, selected)

  curves <- ggplot(
    curves_long,
    aes(x = total, y = probability, color = measure)
  ) +
    geom_hline(
      yintercept = assurance_target,
      linetype = 2,
      linewidth = 0.55
    ) +
    geom_line(linewidth = 0.7) +
    geom_vline(
      data = selected,
      aes(xintercept = total),
      inherit.aes = FALSE,
      color = "#222222",
      linewidth = 0.45
    ) +
    facet_wrap(~round, ncol = 3) +
    scale_color_manual(
      values = c(
        "Raw assurance" = "#176B87",
        "Sustained assurance" = "#D95F02",
        "Sustained simultaneous lower bound" = "#5A3D8A"
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Candidate total sample size",
      y = "Predictive probability",
      title = "Assurance curves from one sequential replay",
      subtitle = "Vertical lines show selected totals"
    )

  calibration <- cluster_bootstrap_calibration(adaptive)

  calibration_plot <- ggplot(
    calibration,
    aes(x = mean_validated_assurance, y = realized_success)
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "#555555"
    ) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.004) +
    geom_point(color = "#176B87", size = 3) +
    coord_equal(xlim = c(0.85, 1.005), ylim = c(0.78, 1.005)) +
    labs(
      x = "Validated assurance",
      y = "Realized width success",
      title = "Assurance calibration"
    )

  expected_realized <- ggplot(
    adaptive,
    aes(
      x = validated_expected_width,
      y = posterior_width,
      color = high_volatility
    )
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "#555555"
    ) +
    geom_point(alpha = 0.55, size = 1.7) +
    scale_color_manual(
      values = c(`FALSE` = "#176B87", `TRUE` = "#D95F02"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    coord_equal() +
    labs(
      x = "Predicted expected width",
      y = "Realized posterior width",
      title = "Predicted and realized precision"
    )

  figure <- curves / (calibration_plot | expected_realized)
  figure <- figure +
    plot_layout(heights = c(1.3, 1)) +
    plot_annotation(title = "Behaviour of the assurance criterion") &
    base_theme
  save_figure(figure, "assurance_behaviour", output_dir, 11, 9)
  calibration
}

plot_filter_performance <- function(
  round_results,
  credibility,
  width_target,
  output_dir
) {
  adaptive <- round_results |>
    filter(policy == "adaptive", feasible) |>
    mutate(error = posterior_mean - true_estimand)

  tracking <- adaptive |>
    group_by(round) |>
    summarise(
      mean_error = mean(error),
      error_q05 = quantile(error, 0.05),
      error_q95 = quantile(error, 0.95),
      coverage = mean(covered),
      width_success = mean(width_success),
      mean_width = mean(posterior_width),
      .groups = "drop"
    )

  error_plot <- ggplot(
    tracking,
    aes(x = round, y = mean_error)
  ) +
    volatility_band() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_ribbon(
      aes(ymin = error_q05, ymax = error_q95),
      fill = "#8DB7C7",
      alpha = 0.45
    ) +
    geom_line(color = "#176B87", linewidth = 0.85) +
    geom_point(color = "#176B87", size = 1.7) +
    labs(
      x = NULL,
      y = "Posterior mean error",
      title = "Filtering error across replicated trajectories",
      subtitle = "Ribbon spans the fifth to ninety-fifth percentiles"
    )

  rates <- tracking |>
    select(round, coverage, width_success) |>
    pivot_longer(
      cols = c(coverage, width_success),
      names_to = "measure",
      values_to = "rate"
    ) |>
    mutate(
      measure = recode(
        measure,
        coverage = "Interval coverage",
        width_success = "Width success"
      )
    )

  rates_plot <- ggplot(
    rates,
    aes(x = round, y = rate, color = measure)
  ) +
    volatility_band() +
    geom_hline(
      yintercept = credibility,
      linetype = 2,
      linewidth = 0.6
    ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.7) +
    scale_color_manual(
      values = c(
        "Interval coverage" = "#8C3B2A",
        "Width success" = "#176B87"
      )
    ) +
    coord_cartesian(ylim = c(0.55, 1.01)) +
    labs(
      x = "Sampling round",
      y = "Empirical rate",
      title = "Coverage and precision attainment",
      subtitle = "Dashed line shows the nominal 0.90 level"
    )

  sample_heatmap <- adaptive |>
    count(round, total_n, name = "replicates")
  heatmap <- ggplot(
    sample_heatmap,
    aes(x = round, y = factor(total_n), fill = replicates)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_viridis_c(option = "C") +
    labs(
      x = "Sampling round",
      y = "Selected total sample size",
      fill = "Replicates",
      title = "Distribution of selected sample sizes"
    ) +
    theme(legend.position = "right")

  precision_plot <- ggplot(
    adaptive,
    aes(
      x = total_n,
      y = posterior_width,
      color = high_volatility
    )
  ) +
    geom_hline(
      yintercept = width_target,
      linetype = 2,
      linewidth = 0.6
    ) +
    geom_point(alpha = 0.55, size = 1.7) +
    scale_color_manual(
      values = c(`FALSE` = "#176B87", `TRUE` = "#D95F02"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    labs(
      x = "Selected total sample size",
      y = "Realized posterior width",
      title = "Effort and realized precision"
    )

  figure <- (error_plot / rates_plot) | (heatmap / precision_plot)
  figure <- figure +
    plot_annotation(title = "Sequential filtering and design performance") &
    base_theme
  save_figure(figure, "filter_and_design_performance", output_dir, 11, 8)
  tracking
}

plot_policy_comparison <- function(round_results, output_dir) {
  replicate_summary <- round_results |>
    filter(feasible) |>
    group_by(replicate, policy) |>
    summarise(
      total_effort = sum(total_n),
      mean_width = mean(posterior_width),
      width_success = mean(width_success),
      coverage = mean(covered),
      rmse = sqrt(mean(squared_error)),
      .groups = "drop"
    )

  effort <- ggplot(
    replicate_summary,
    aes(x = policy, y = total_effort, fill = policy)
  ) +
    geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.75) +
    geom_point(
      position = position_jitter(width = 0.08),
      alpha = 0.65,
      size = 1.5
    ) +
    scale_fill_manual(
      values = c(adaptive = "#176B87", fixed_total = "#8C3B2A"),
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    scale_x_discrete(
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    labs(
      x = NULL,
      y = "Total samples over 12 rounds",
      title = "Sampling effort by replicate"
    ) +
    guides(fill = "none")

  quality <- replicate_summary |>
    select(replicate, policy, width_success, coverage) |>
    pivot_longer(
      cols = c(width_success, coverage),
      names_to = "measure",
      values_to = "rate"
    ) |>
    mutate(
      measure = recode(
        measure,
        width_success = "Width success",
        coverage = "Interval coverage"
      )
    )

  quality_plot <- ggplot(
    quality,
    aes(x = policy, y = rate, color = policy)
  ) +
    geom_hline(yintercept = 0.90, linetype = 2, linewidth = 0.6) +
    geom_boxplot(width = 0.48, outlier.shape = NA) +
    geom_point(
      position = position_jitter(width = 0.08),
      alpha = 0.65,
      size = 1.5
    ) +
    facet_wrap(~measure) +
    scale_color_manual(
      values = c(adaptive = "#176B87", fixed_total = "#8C3B2A"),
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    scale_x_discrete(
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    coord_cartesian(ylim = c(0.55, 1.01)) +
    labs(
      x = NULL,
      y = "Replicate-level rate",
      title = "Precision and coverage by replicate"
    ) +
    guides(color = "none")

  efficiency <- round_results |>
    filter(feasible) |>
    ggplot(
      aes(
        x = total_n,
        y = posterior_width,
        color = policy
      )
    ) +
    geom_point(alpha = 0.38, size = 1.5) +
    scale_color_manual(
      values = c(adaptive = "#176B87", fixed_total = "#8C3B2A"),
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    labs(
      x = "Round sample size",
      y = "Posterior interval width",
      title = "Round-level effort and precision"
    )

  figure <- effort | quality_plot / efficiency
  figure <- figure +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(title = "Adaptive and fixed-total policy comparison") &
    base_theme
  save_figure(figure, "policy_comparison", output_dir, 11, 8)
  replicate_summary
}

plot_aggregate_assurance <- function(
  assurance_curves,
  assurance_target,
  output_dir
) {
  curve_summary <- assurance_curves |>
    group_by(round, high_volatility, total) |>
    summarise(
      raw_mean = mean(assurance),
      raw_q10 = quantile(assurance, 0.10),
      raw_q90 = quantile(assurance, 0.90),
      sustained_mean = mean(sustained_assurance),
      sustained_q10 = quantile(sustained_assurance, 0.10),
      sustained_q90 = quantile(sustained_assurance, 0.90),
      sustained_lower_mean = mean(sustained_lower_bound),
      sustained_lower_q10 = quantile(sustained_lower_bound, 0.10),
      sustained_lower_q90 = quantile(sustained_lower_bound, 0.90),
      selection_rate = mean(selected),
      mean_n_A = mean(n_A),
      mean_n_B = mean(n_B),
      .groups = "drop"
    )

  display_rounds <- c(1, 5, 6, 7, 8, 12)
  displayed <- filter(curve_summary, round %in% display_rounds)
  curves <- ggplot(displayed, aes(x = total)) +
    geom_hline(
      yintercept = assurance_target,
      linetype = 2,
      linewidth = 0.55
    ) +
    geom_ribbon(
      aes(ymin = raw_q10, ymax = raw_q90),
      fill = "#8DB7C7",
      alpha = 0.32
    ) +
    geom_line(
      aes(y = raw_mean, color = "Raw assurance"),
      linewidth = 0.75
    ) +
    geom_line(
      aes(y = sustained_mean, color = "Sustained assurance"),
      linewidth = 0.75,
      linetype = 3
    ) +
    geom_line(
      aes(
        y = sustained_lower_mean,
        color = "Sustained simultaneous lower bound"
      ),
      linewidth = 0.75,
      linetype = 4
    ) +
    facet_wrap(~round, ncol = 3) +
    scale_color_manual(
      values = c(
        "Raw assurance" = "#176B87",
        "Sustained assurance" = "#D95F02",
        "Sustained simultaneous lower bound" = "#5A3D8A"
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Candidate total sample size",
      y = "Assurance",
      title = "Assurance across replicated filtering histories",
      subtitle = "Ribbons span the tenth to ninetieth percentiles"
    )

  selected <- assurance_curves |>
    filter(selected)
  selected_sizes <- ggplot(
    selected,
    aes(
      x = factor(round),
      y = total,
      fill = high_volatility
    )
  ) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.30) +
    scale_fill_manual(
      values = c(`FALSE` = "#8DB7C7", `TRUE` = "#E9A36A"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    labs(
      x = "Sampling round",
      y = "Selected total sample size",
      title = "Distribution of selected designs"
    )

  curve_diagnostics <- assurance_curves |>
    arrange(replicate, round, total) |>
    group_by(replicate, round, high_volatility) |>
    summarise(
      nonmonotone = any(diff(assurance) < 0),
      maximum_downward_step = max(c(0, -diff(assurance))),
      maximum_sustained_gap =
        max(assurance - sustained_assurance),
      .groups = "drop"
    )

  nonmonotonicity <- curve_diagnostics |>
    group_by(round, high_volatility) |>
    summarise(
      proportion_nonmonotone = mean(nonmonotone),
      mean_maximum_gap = mean(maximum_sustained_gap),
      .groups = "drop"
    ) |>
    ggplot(
      aes(
        x = factor(round),
        y = proportion_nonmonotone,
        fill = high_volatility
      )
    ) +
    geom_col(width = 0.65) +
    scale_fill_manual(
      values = c(`FALSE` = "#8DB7C7", `TRUE` = "#E9A36A"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Sampling round",
      y = "Proportion of curves",
      title = "Observed nonmonotonic assurance",
      subtitle = "Candidate totals share predictive Monte Carlo streams"
    )

  figure <- curves / (selected_sizes | nonmonotonicity)
  figure <- figure +
    plot_layout(heights = c(1.35, 1)) +
    plot_annotation(title = "Aggregate assurance behaviour") &
    base_theme
  save_figure(figure, "aggregate_assurance_behaviour", output_dir, 11, 9)

  list(
    curve_summary = curve_summary,
    curve_diagnostics = curve_diagnostics
  )
}

plot_numerical_diagnostics <- function(
  round_results,
  assurance_target,
  output_dir
) {
  adaptive <- round_results |>
    filter(policy == "adaptive", feasible) |>
    mutate(
      assurance_difference =
        predicted_assurance - validated_assurance
    )

  validation <- ggplot(
    adaptive,
    aes(
      x = predicted_assurance,
      y = validated_assurance,
      color = high_volatility
    )
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "#555555"
    ) +
    geom_hline(
      yintercept = assurance_target,
      linetype = 3,
      linewidth = 0.5
    ) +
    geom_point(alpha = 0.45, size = 1.5) +
    scale_color_manual(
      values = c(`FALSE` = "#176B87", `TRUE` = "#D95F02"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    coord_equal(xlim = c(0.88, 1.005), ylim = c(0.88, 1.005)) +
    labs(
      x = "Design-stage assurance estimate",
      y = "Independent validation estimate",
      title = "Selected-design validation"
    )

  difference <- ggplot(
    adaptive,
    aes(x = assurance_difference, fill = high_volatility)
  ) +
    geom_histogram(
      bins = 35,
      position = "identity",
      alpha = 0.58
    ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6) +
    scale_fill_manual(
      values = c(`FALSE` = "#8DB7C7", `TRUE` = "#E9A36A"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    labs(
      x = "Design estimate minus validation estimate",
      y = "Decisions",
      title = "Design and validation differences"
    )

  ess <- ggplot(
    adaptive,
    aes(
      x = factor(round),
      y = ess_fraction,
      fill = resampled
    )
  ) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.25) +
    geom_hline(yintercept = 0.50, linetype = 2, linewidth = 0.6) +
    scale_fill_manual(
      values = c(`FALSE` = "#8DB7C7", `TRUE` = "#E9A36A"),
      labels = c(`FALSE` = "Not resampled", `TRUE` = "Resampled")
    ) +
    labs(
      x = "Sampling round",
      y = "Effective sample size fraction",
      title = "Particle-filter stability",
      subtitle = "Dashed line shows the resampling threshold"
    )

  mc_error <- ggplot(
    adaptive,
    aes(
      x = validation_mc_se,
      y = abs(assurance_difference),
      color = high_volatility
    )
  ) +
    geom_point(alpha = 0.45, size = 1.5) +
    scale_color_manual(
      values = c(`FALSE` = "#176B87", `TRUE` = "#D95F02"),
      labels = c(`FALSE` = "Ordinary rounds", `TRUE` = "High volatility")
    ) +
    labs(
      x = "Validation Monte Carlo standard error",
      y = "Absolute design-validation difference",
      title = "Monte Carlo uncertainty"
    )

  figure <- (validation | difference) / (ess | mc_error)
  figure <- figure +
    plot_annotation(title = "Numerical diagnostics") &
    base_theme
  save_figure(figure, "numerical_diagnostics", output_dir, 11, 8)

  adaptive |>
    summarise(
      decisions = n(),
      mean_design_assurance = mean(predicted_assurance),
      mean_sustained_lower_bound =
        mean(sustained_assurance_lower_bound),
      mean_validation_assurance = mean(validated_assurance),
      mean_design_validation_difference =
        mean(assurance_difference),
      mean_absolute_design_validation_difference =
        mean(abs(assurance_difference)),
      rmse_design_validation_difference =
        sqrt(mean(assurance_difference^2)),
      selected_below_target_after_validation =
        mean(validated_assurance < assurance_target),
      mean_design_mc_se = mean(assurance_mc_se),
      mean_validation_mc_se = mean(validation_mc_se),
      mean_ess_fraction = mean(ess_fraction),
      resampling_rate = mean(resampled)
    )
}

summarise_framework_components <- function(
  round_results,
  curve_diagnostics,
  assurance_target
) {
  adaptive <- round_results |>
    filter(policy == "adaptive", feasible) |>
    mutate(
      error = posterior_mean - true_estimand,
      allocated_share_A = n_A / total_n,
      proxy_share_A =
        sqrt(allocation_coefficient_A) /
        (
          sqrt(allocation_coefficient_A) +
            sqrt(allocation_coefficient_B)
        )
    )
  fixed <- round_results |>
    filter(policy == "fixed_total", feasible)

  adaptive_rep <- adaptive |>
    group_by(replicate) |>
    summarise(
      predictive_coverage = mean(predictive_covered),
      posterior_coverage = mean(covered),
      rmse = sqrt(mean(squared_error)),
      width_success = mean(width_success),
      selected_below_target =
        mean(.data$validated_assurance < assurance_target),
      mean_validated_assurance = mean(.data$validated_assurance),
      mean_sustained_lower_bound =
        mean(sustained_assurance_lower_bound),
      validation_calibration =
        width_success - mean_validated_assurance,
      ordinary_mean_n = mean(total_n[!high_volatility]),
      volatile_mean_n = mean(total_n[high_volatility]),
      uncertainty_effort_correlation =
        cor(predictive_width, total_n, method = "spearman"),
      allocation_share_mad =
        mean(abs(allocated_share_A - proxy_share_A)),
      mean_ess_fraction = mean(ess_fraction),
      resampling_rate = mean(resampled),
      design_validation_mae =
        mean(abs(
          predicted_assurance - .data$validated_assurance
        )),
      total_effort = sum(total_n),
      .groups = "drop"
    )
  fixed_rep <- fixed |>
    group_by(replicate) |>
    summarise(
      fixed_total_effort = sum(total_n),
      fixed_width_success = mean(width_success),
      fixed_coverage = mean(covered),
      .groups = "drop"
    )
  curve_rep <- curve_diagnostics |>
    group_by(replicate) |>
    summarise(
      nonmonotone_round_fraction = mean(nonmonotone),
      mean_maximum_sustained_gap = mean(maximum_sustained_gap),
      .groups = "drop"
    )

  replicate_metrics <- adaptive_rep |>
    left_join(fixed_rep, by = "replicate") |>
    left_join(curve_rep, by = "replicate") |>
    mutate(
      paired_effort_reduction =
        1 - total_effort / fixed_total_effort
    )

  definitions <- data.frame(
    component = c(
      "Filtering",
      "Filtering",
      "Filtering",
      "Design response",
      "Design response",
      "Design response",
      "Allocation",
      "Assurance",
      "Assurance",
      "Assurance",
      "Assurance",
      "Assurance",
      "Sustained assurance",
      "Sustained assurance",
      "Numerical",
      "Numerical",
      "Efficiency"
    ),
    metric = c(
      "Predictive interval coverage",
      "Posterior interval coverage",
      "Posterior RMSE",
      "Mean sample size in ordinary rounds",
      "Mean sample size in high-volatility rounds",
      "Predictive-width and effort correlation",
      "Allocation-share mean absolute deviation",
      "Validated assurance",
      "Sustained simultaneous lower bound",
      "Realized width success",
      "Validation calibration difference",
      "Selected designs below target after validation",
      "Fraction of nonmonotone round curves",
      "Mean maximum raw-sustained gap",
      "Mean ESS fraction",
      "Design-validation assurance error",
      "Total effort over 12 rounds"
    ),
    variable = c(
      "predictive_coverage",
      "posterior_coverage",
      "rmse",
      "ordinary_mean_n",
      "volatile_mean_n",
      "uncertainty_effort_correlation",
      "allocation_share_mad",
      "mean_validated_assurance",
      "mean_sustained_lower_bound",
      "width_success",
      "validation_calibration",
      "selected_below_target",
      "nonmonotone_round_fraction",
      "mean_maximum_sustained_gap",
      "mean_ess_fraction",
      "design_validation_mae",
      "total_effort"
    ),
    stringsAsFactors = FALSE
  )

  summaries <- lapply(seq_len(nrow(definitions)), function(index) {
    variable <- definitions$variable[[index]]
    values <- replicate_metrics[[variable]]
    standard_error <- sd(values) / sqrt(length(values))
    lower <- mean(values) - 1.96 * standard_error
    upper <- mean(values) + 1.96 * standard_error
    probability_metrics <- c(
      "predictive_coverage",
      "posterior_coverage",
      "width_success",
      "mean_validated_assurance",
      "mean_sustained_lower_bound",
      "selected_below_target",
      "nonmonotone_round_fraction",
      "mean_ess_fraction",
      "resampling_rate"
    )
    nonnegative_metrics <- c(
      "rmse",
      "ordinary_mean_n",
      "volatile_mean_n",
      "allocation_share_mad",
      "mean_maximum_sustained_gap",
      "design_validation_mae",
      "total_effort"
    )
    if (variable %in% probability_metrics) {
      lower <- max(0, lower)
      upper <- min(1, upper)
    } else if (variable == "uncertainty_effort_correlation") {
      lower <- max(-1, lower)
      upper <- min(1, upper)
    } else if (variable %in% nonnegative_metrics) {
      lower <- max(0, lower)
    }
    data.frame(
      component = definitions$component[[index]],
      metric = definitions$metric[[index]],
      estimate = mean(values),
      standard_error = standard_error,
      lower_95 = lower,
      upper_95 = upper,
      replicates = length(values)
    )
  })
  list(
    summary = bind_rows(summaries),
    replicate_metrics = replicate_metrics
  )
}

write_framework_latex_table <- function(summary, output_dir) {
  escape_latex <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([&_#%])", "\\\\\\1", x)
    x
  }
  format_value <- function(metric, value) {
    count_metrics <- c(
      "Mean sample size in ordinary rounds",
      "Mean sample size in high-volatility rounds",
      "Total effort over 12 rounds"
    )
    if (metric %in% count_metrics) {
      sprintf("%.1f", value)
    } else {
      sprintf("%.3f", value)
    }
  }

  rows <- vapply(seq_len(nrow(summary)), function(index) {
    metric <- summary$metric[[index]]
    estimate <- format_value(metric, summary$estimate[[index]])
    interval <- paste0(
      "(",
      format_value(metric, summary$lower_95[[index]]),
      ", ",
      format_value(metric, summary$upper_95[[index]]),
      ")"
    )
    paste(
      escape_latex(summary$component[[index]]),
      escape_latex(metric),
      estimate,
      interval,
      sep = " & "
    ) |>
      paste0(" \\\\")
  }, character(1))

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Simulation summary across replicated trajectories. Intervals are normal-approximation 95\\% Monte Carlo intervals for the replicate-level mean.}",
    "\\label{tab:simple_simulation_summary}",
    "\\footnotesize",
    "\\begin{tabular}{llrr}",
    "\\toprule",
    "Component & Measure & Estimate & 95\\% interval \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(
    lines,
    file.path(output_dir, "framework_component_summary_table.tex")
  )
}

main <- function() {
  input_dir <- parse_input_dir(commandArgs(trailingOnly = TRUE))
  rounds <- read.csv(
    file.path(input_dir, "simulation_rounds.csv"),
    stringsAsFactors = FALSE
  )
  assurance_curves <- read.csv(
    file.path(input_dir, "assurance_curves.csv"),
    stringsAsFactors = FALSE
  )
  config <- read.csv(
    file.path(input_dir, "simulation_config.csv"),
    stringsAsFactors = FALSE
  )

  width_target <- read_config_value(config, "width target")
  assurance_target <- read_config_value(config, "assurance target")
  credibility <- read_config_value(config, "credible level")
  adaptive <- rounds |>
    filter(policy == "adaptive", feasible)

  dynamic_summary <- plot_dynamic_response(
    adaptive,
    width_target,
    input_dir
  )
  allocation_summary <- plot_allocation_behaviour(adaptive, input_dir)
  calibration_summary <- plot_assurance_behaviour(
    adaptive,
    assurance_curves,
    assurance_target,
    input_dir
  )
  tracking_summary <- plot_filter_performance(
    rounds,
    credibility,
    width_target,
    input_dir
  )
  replicate_summary <- plot_policy_comparison(rounds, input_dir)
  aggregate_assurance <- plot_aggregate_assurance(
    assurance_curves,
    assurance_target,
    input_dir
  )
  numerical_summary <- plot_numerical_diagnostics(
    rounds,
    assurance_target,
    input_dir
  )
  framework_summary <- summarise_framework_components(
    rounds,
    aggregate_assurance$curve_diagnostics,
    assurance_target
  )

  write.csv(
    dynamic_summary,
    file.path(input_dir, "dynamic_round_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    allocation_summary,
    file.path(input_dir, "allocation_round_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    calibration_summary,
    file.path(input_dir, "assurance_calibration_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    tracking_summary,
    file.path(input_dir, "filter_tracking_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    replicate_summary,
    file.path(input_dir, "replicate_policy_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    aggregate_assurance$curve_summary,
    file.path(input_dir, "candidate_curve_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    aggregate_assurance$curve_diagnostics,
    file.path(input_dir, "curve_diagnostics.csv"),
    row.names = FALSE
  )
  write.csv(
    numerical_summary,
    file.path(input_dir, "numerical_diagnostics_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    framework_summary$summary,
    file.path(input_dir, "framework_component_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    framework_summary$replicate_metrics,
    file.path(input_dir, "framework_replicate_metrics.csv"),
    row.names = FALSE
  )
  write_framework_latex_table(framework_summary$summary, input_dir)

  message("Generated behaviour figures in ", input_dir)
}

main()
