suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) plogis(x)

parse_options <- function(args) {
  defaults <- list(
    replicates = 20L,
    particles = 400L,
    design_mc = 250L,
    validation_mc = 1000L,
    workers = 1L,
    selection_confidence = 0.95,
    seed = 20260729L,
    output_dir = file.path(
      "output",
      "simple_framework_simulation"
    )
  )

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Options must use --name=value syntax")
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    name <- gsub("-", "_", pieces[[1]])
    value <- paste(pieces[-1], collapse = "=")
    if (!name %in% names(defaults)) {
      stop("Unknown option --", pieces[[1]])
    }
    defaults[[name]] <- if (is.character(defaults[[name]])) {
      value
    } else if (is.integer(defaults[[name]])) {
      as.integer(value)
    } else {
      as.numeric(value)
    }
  }
  defaults
}

rmvn <- function(n, mean, covariance) {
  z <- matrix(rnorm(n * length(mean)), nrow = n)
  sweep(z %*% chol(covariance), 2, mean, FUN = "+")
}

normalise_log_weights <- function(log_weights) {
  maximum <- max(log_weights)
  weights <- exp(log_weights - maximum)
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("Particle weights could not be normalised")
  }
  weights / total
}

systematic_resample <- function(weights) {
  n <- length(weights)
  positions <- (runif(1) + 0:(n - 1)) / n
  cumulative <- pmin(cumsum(weights), 1)
  cumulative[[n]] <- 1
  findInterval(positions, c(0, cumulative), rightmost.closed = TRUE)
}

weighted_interval <- function(values, weights, credibility) {
  ordering <- order(values)
  values <- values[ordering]
  cumulative <- cumsum(weights[ordering]) / sum(weights)
  tail_probability <- (1 - credibility) / 2
  lower <- values[which(cumulative >= tail_probability - 1e-12)[1]]
  upper <- values[which(cumulative >= 1 - tail_probability - 1e-12)[1]]
  c(lower = lower, upper = upper, width = upper - lower)
}

particle_log_likelihood <- function(probabilities, counts, allocation) {
  probabilities <- pmin(pmax(probabilities, 1e-10), 1 - 1e-10)
  result <- rep(0, nrow(probabilities))
  for (s in seq_along(allocation)) {
    result <- result +
      counts[[s]] * log(probabilities[, s]) +
      (allocation[[s]] - counts[[s]]) * log1p(-probabilities[, s])
  }
  result
}

simultaneous_binomial_lower_bound <- function(
  successes,
  trials,
  family_confidence,
  comparisons
) {
  if (trials <= 0L || comparisons <= 0L) {
    stop("Trials and comparisons must be positive")
  }
  if (family_confidence <= 0 || family_confidence >= 1) {
    stop("Selection confidence must lie strictly between zero and one")
  }
  tail_probability <- (1 - family_confidence) / comparisons
  if (successes <= 0L) {
    return(0)
  }
  qbeta(
    tail_probability,
    shape1 = successes,
    shape2 = trials - successes + 1
  )
}

make_config <- function(options) {
  initial_sd <- c(0.30, 0.30)
  process_sd <- c(0.10, 0.13)
  correlation <- 0.45
  correlation_matrix <- matrix(
    c(1, correlation, correlation, 1),
    nrow = 2
  )

  list(
    strata = c("A", "B"),
    stratum_weights = c(0.60, 0.40),
    initial_mean = logit(c(0.12, 0.28)),
    initial_covariance =
      diag(initial_sd) %*% correlation_matrix %*% diag(initial_sd),
    process_covariance =
      diag(process_sd) %*% correlation_matrix %*% diag(process_sd),
    rounds = 12L,
    volatility_multiplier = c(1, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 1),
    particles = options$particles,
    resampling_fraction = 0.50,
    credibility = 0.90,
    width_target = 0.10,
    assurance_target = 0.90,
    stratum_minimum = c(5L, 5L),
    candidate_totals = seq.int(20L, 240L, by = 20L),
    fixed_total = 240L,
    design_mc = options$design_mc,
    validation_mc = options$validation_mc,
    selection_confidence = options$selection_confidence,
    replicates = options$replicates,
    workers = options$workers,
    seed = options$seed,
    output_dir = options$output_dir
  )
}

simulate_truth <- function(config) {
  states <- matrix(
    NA_real_,
    nrow = config$rounds + 1L,
    ncol = length(config$strata),
    dimnames = list(0:config$rounds, config$strata)
  )
  states[1, ] <- rmvn(1, config$initial_mean, config$initial_covariance)
  for (t in seq_len(config$rounds)) {
    covariance <- config$process_covariance *
      config$volatility_multiplier[[t]]^2
    states[t + 1L, ] <- states[t, ] + rmvn(
      1,
      rep(0, length(config$strata)),
      covariance
    )
  }
  states
}

propagate_particles <- function(states, round_index, config) {
  covariance <- config$process_covariance *
    config$volatility_multiplier[[round_index]]^2
  states + rmvn(
    nrow(states),
    rep(0, length(config$strata)),
    covariance
  )
}

build_greedy_family <- function(
  probabilities,
  predictive_weights,
  config
) {
  coefficients <- config$stratum_weights^2 *
    colSums(
      probabilities * (1 - probabilities) * predictive_weights
    )
  allocation <- as.integer(config$stratum_minimum)
  maximum_total <- max(config$candidate_totals)
  allocations <- matrix(
    NA_integer_,
    nrow = length(config$candidate_totals),
    ncol = length(config$strata),
    dimnames = list(config$candidate_totals, config$strata)
  )

  for (total in seq.int(sum(allocation), maximum_total)) {
    if (total %in% config$candidate_totals) {
      allocations[as.character(total), ] <- allocation
    }
    if (total < maximum_total) {
      marginal_gain <- coefficients / (allocation * (allocation + 1))
      chosen_stratum <- which.max(marginal_gain)
      allocation[[chosen_stratum]] <- allocation[[chosen_stratum]] + 1L
    }
  }

  if (anyNA(allocations) ||
      any(rowSums(allocations) != config$candidate_totals)) {
    stop("The greedy allocation family is internally inconsistent")
  }
  allocations
}

hypothetical_widths <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation,
  prospective_counts,
  credibility
) {
  n_particles <- nrow(probabilities)
  n_draws <- nrow(prospective_counts)
  log_probabilities <- log(pmin(pmax(probabilities, 1e-10), 1 - 1e-10))
  log_complements <- log1p(
    -pmin(pmax(probabilities, 1e-10), 1 - 1e-10)
  )

  log_weight_matrix <- matrix(
    log(pmax(predictive_weights, 1e-300)),
    nrow = n_particles,
    ncol = n_draws
  )
  for (s in seq_along(allocation)) {
    log_weight_matrix <- log_weight_matrix +
      tcrossprod(log_probabilities[, s], prospective_counts[, s]) +
      tcrossprod(
        log_complements[, s],
        allocation[[s]] - prospective_counts[, s]
      )
  }

  column_maxima <- apply(log_weight_matrix, 2, max)
  posterior_weights <- exp(
    sweep(log_weight_matrix, 2, column_maxima, FUN = "-")
  )
  posterior_weights <- sweep(
    posterior_weights,
    2,
    colSums(posterior_weights),
    FUN = "/"
  )

  ordering <- order(estimands)
  sorted_estimands <- estimands[ordering]
  cumulative <- apply(posterior_weights[ordering, , drop = FALSE], 2, cumsum)
  tail_probability <- (1 - credibility) / 2
  lower_index <- max.col(
    t(cumulative >= tail_probability),
    ties.method = "first"
  )
  upper_index <- max.col(
    t(cumulative >= 1 - tail_probability),
    ties.method = "first"
  )
  sorted_estimands[upper_index] - sorted_estimands[lower_index]
}

make_predictive_mc_stream <- function(
  probabilities,
  predictive_weights,
  maximum_allocation,
  mc_draws
) {
  particle_index <- sample.int(
    nrow(probabilities),
    mc_draws,
    replace = TRUE,
    prob = predictive_weights
  )
  cumulative_counts <- lapply(
    seq_along(maximum_allocation),
    function(s) {
      maximum <- maximum_allocation[[s]]
      if (maximum <= 0L) {
        return(matrix(0L, nrow = mc_draws, ncol = 0L))
      }
      uniforms <- matrix(
        runif(mc_draws * maximum),
        nrow = mc_draws,
        ncol = maximum
      )
      successes <- uniforms <= probabilities[particle_index, s]
      cumulative <- matrix(
        0L,
        nrow = mc_draws,
        ncol = maximum
      )
      cumulative[, 1] <- successes[, 1]
      if (maximum > 1L) {
        for (sample_index in 2:maximum) {
          cumulative[, sample_index] <-
            cumulative[, sample_index - 1L] +
            successes[, sample_index]
        }
      }
      cumulative
    }
  )
  list(
    particle_index = particle_index,
    cumulative_counts = cumulative_counts,
    mc_draws = mc_draws
  )
}

counts_from_mc_stream <- function(stream, allocation) {
  counts <- matrix(
    0L,
    nrow = stream$mc_draws,
    ncol = length(allocation)
  )
  for (s in seq_along(allocation)) {
    if (allocation[[s]] > 0L) {
      counts[, s] <-
        stream$cumulative_counts[[s]][, allocation[[s]]]
    }
  }
  counts
}

estimate_assurance <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation,
  config,
  mc_draws = config$design_mc,
  prospective_counts = NULL,
  comparisons = 1L
) {
  if (is.null(prospective_counts)) {
    stream <- make_predictive_mc_stream(
      probabilities,
      predictive_weights,
      allocation,
      mc_draws
    )
    prospective_counts <- counts_from_mc_stream(stream, allocation)
  } else {
    if (nrow(prospective_counts) != mc_draws ||
        ncol(prospective_counts) != length(allocation)) {
      stop("Prospective counts do not match the assurance calculation")
    }
  }
  widths <- hypothetical_widths(
    probabilities,
    estimands,
    predictive_weights,
    allocation,
    prospective_counts,
    config$credibility
  )
  successes <- sum(widths <= config$width_target)
  assurance <- successes / mc_draws
  c(
    assurance = assurance,
    assurance_lower_bound = simultaneous_binomial_lower_bound(
      successes,
      mc_draws,
      config$selection_confidence,
      comparisons
    ),
    mc_se = sqrt(assurance * (1 - assurance) / mc_draws),
    expected_width = mean(widths)
  )
}

select_adaptive_design <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation_family,
  config
) {
  candidate_count <- nrow(allocation_family)
  design_stream <- make_predictive_mc_stream(
    probabilities,
    predictive_weights,
    apply(allocation_family, 2, max),
    config$design_mc
  )
  evaluations <- t(vapply(
    seq_len(candidate_count),
    function(index) {
      estimate_assurance(
        probabilities,
        estimands,
        predictive_weights,
        allocation_family[index, ],
        config,
        prospective_counts = counts_from_mc_stream(
          design_stream,
          allocation_family[index, ]
        ),
        comparisons = candidate_count
      )
    },
    numeric(4)
  ))
  raw_assurance <- evaluations[, "assurance"]
  assurance_lower_bound <- evaluations[, "assurance_lower_bound"]
  sustained_assurance <- rev(cummin(rev(raw_assurance)))
  sustained_lower_bound <- rev(cummin(rev(assurance_lower_bound)))
  curve <- cbind(
    total = config$candidate_totals,
    evaluations,
    sustained_assurance = sustained_assurance,
    sustained_lower_bound = sustained_lower_bound
  )
  feasible <- which(sustained_lower_bound >= config$assurance_target)
  if (length(feasible) == 0L) {
    return(list(feasible = FALSE, curve = curve))
  }

  selected <- feasible[[1]]
  list(
    feasible = TRUE,
    total = config$candidate_totals[[selected]],
    allocation = allocation_family[selected, ],
    raw_assurance = raw_assurance[[selected]],
    sustained_assurance = sustained_assurance[[selected]],
    assurance_lower_bound = assurance_lower_bound[[selected]],
    sustained_lower_bound = sustained_lower_bound[[selected]],
    mc_se = evaluations[selected, "mc_se"],
    expected_width = evaluations[selected, "expected_width"],
    curve = curve
  )
}

select_fixed_total_design <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation_family,
  config
) {
  selected <- match(config$fixed_total, config$candidate_totals)
  allocation <- allocation_family[selected, ]
  evaluation <- estimate_assurance(
    probabilities,
    estimands,
    predictive_weights,
    allocation,
    config
  )
  list(
    feasible = TRUE,
    total = config$fixed_total,
    allocation = allocation,
    raw_assurance = unname(evaluation[["assurance"]]),
    sustained_assurance = NA_real_,
    assurance_lower_bound =
      unname(evaluation[["assurance_lower_bound"]]),
    sustained_lower_bound = NA_real_,
    mc_se = unname(evaluation[["mc_se"]]),
    expected_width = unname(evaluation[["expected_width"]])
  )
}

run_policy <- function(truth_states, policy, replicate_id, config) {
  policy_offset <- if (policy == "adaptive") 100000L else 200000L
  set.seed(config$seed + replicate_id * 1000L + policy_offset)

  particles <- rmvn(
    config$particles,
    config$initial_mean,
    config$initial_covariance
  )
  weights <- rep(1 / config$particles, config$particles)
  rows <- vector("list", config$rounds)
  curve_rows <- list()

  for (t in seq_len(config$rounds)) {
    particles <- propagate_particles(particles, t, config)
    probabilities <- inv_logit(particles)
    estimands <- as.vector(probabilities %*% config$stratum_weights)
    predictive_interval <- weighted_interval(
      estimands,
      weights,
      config$credibility
    )
    predictive_mean <- sum(weights * estimands)
    predictive_stratum_means <- colSums(probabilities * weights)
    allocation_coefficients <- config$stratum_weights^2 *
      colSums(probabilities * (1 - probabilities) * weights)
    allocation_family <- build_greedy_family(
      probabilities,
      weights,
      config
    )

    design <- if (policy == "adaptive") {
      select_adaptive_design(
        probabilities,
        estimands,
        weights,
        allocation_family,
        config
      )
    } else {
      select_fixed_total_design(
        probabilities,
        estimands,
        weights,
        allocation_family,
        config
      )
    }

    if (policy == "adaptive") {
      curve <- as.data.frame(design$curve)
      curve$replicate <- replicate_id
      curve$round <- t
      curve$high_volatility <- config$volatility_multiplier[[t]] > 1
      curve$predictive_width <- predictive_interval[["width"]]
      curve$n_A <- allocation_family[, 1]
      curve$n_B <- allocation_family[, 2]
      curve$selected <- if (design$feasible) {
        curve$total == design$total
      } else {
        FALSE
      }
      curve_rows[[length(curve_rows) + 1L]] <- curve
    }

    if (!design$feasible) {
      rows[[t]] <- data.frame(
        replicate = replicate_id,
        policy = policy,
        round = t,
        high_volatility = config$volatility_multiplier[[t]] > 1,
        volatility_multiplier = config$volatility_multiplier[[t]],
        feasible = FALSE,
        true_p_A = NA_real_,
        true_p_B = NA_real_,
        true_estimand = NA_real_,
        predictive_p_A = predictive_stratum_means[[1]],
        predictive_p_B = predictive_stratum_means[[2]],
        predictive_mean = predictive_mean,
        predictive_lower = predictive_interval[["lower"]],
        predictive_upper = predictive_interval[["upper"]],
        predictive_width = predictive_interval[["width"]],
        predictive_covered = NA,
        allocation_coefficient_A = allocation_coefficients[[1]],
        allocation_coefficient_B = allocation_coefficients[[2]],
        n_A = NA_integer_,
        n_B = NA_integer_,
        total_n = NA_integer_,
        y_A = NA_integer_,
        y_B = NA_integer_,
        predicted_assurance = NA_real_,
        sustained_assurance = NA_real_,
        assurance_lower_bound = NA_real_,
        sustained_assurance_lower_bound = NA_real_,
        assurance_mc_se = NA_real_,
        predicted_expected_width = NA_real_,
        validated_assurance = NA_real_,
        validation_mc_se = NA_real_,
        validated_expected_width = NA_real_,
        posterior_mean = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        posterior_width = NA_real_,
        width_success = NA,
        covered = NA,
        squared_error = NA_real_,
        effective_sample_size = NA_real_,
        ess_fraction = NA_real_,
        resampled = NA
      )
      break
    }

    validation <- estimate_assurance(
      probabilities,
      estimands,
      weights,
      design$allocation,
      config,
      mc_draws = config$validation_mc
    )
    true_probabilities <- inv_logit(truth_states[t + 1L, ])
    true_estimand <- sum(config$stratum_weights * true_probabilities)
    observed_counts <- vapply(
      seq_along(design$allocation),
      function(s) {
        rbinom(1, design$allocation[[s]], true_probabilities[[s]])
      },
      integer(1)
    )

    log_weights <- log(pmax(weights, 1e-300)) +
      particle_log_likelihood(
        probabilities,
        observed_counts,
        design$allocation
      )
    weights <- normalise_log_weights(log_weights)
    interval <- weighted_interval(
      estimands,
      weights,
      config$credibility
    )
    posterior_mean <- sum(weights * estimands)
    effective_sample_size <- 1 / sum(weights^2)

    rows[[t]] <- data.frame(
      replicate = replicate_id,
      policy = policy,
      round = t,
      high_volatility = config$volatility_multiplier[[t]] > 1,
      volatility_multiplier = config$volatility_multiplier[[t]],
      feasible = TRUE,
      true_p_A = true_probabilities[[1]],
      true_p_B = true_probabilities[[2]],
      true_estimand = true_estimand,
      predictive_p_A = predictive_stratum_means[[1]],
      predictive_p_B = predictive_stratum_means[[2]],
      predictive_mean = predictive_mean,
      predictive_lower = predictive_interval[["lower"]],
      predictive_upper = predictive_interval[["upper"]],
      predictive_width = predictive_interval[["width"]],
      predictive_covered =
        predictive_interval[["lower"]] <= true_estimand &&
          true_estimand <= predictive_interval[["upper"]],
      allocation_coefficient_A = allocation_coefficients[[1]],
      allocation_coefficient_B = allocation_coefficients[[2]],
      n_A = design$allocation[[1]],
      n_B = design$allocation[[2]],
      total_n = design$total,
      y_A = observed_counts[[1]],
      y_B = observed_counts[[2]],
      predicted_assurance = design$raw_assurance,
      sustained_assurance = design$sustained_assurance,
      assurance_lower_bound = design$assurance_lower_bound,
      sustained_assurance_lower_bound = design$sustained_lower_bound,
      assurance_mc_se = design$mc_se,
      predicted_expected_width = design$expected_width,
      validated_assurance = unname(validation[["assurance"]]),
      validation_mc_se = unname(validation[["mc_se"]]),
      validated_expected_width =
        unname(validation[["expected_width"]]),
      posterior_mean = posterior_mean,
      lower = interval[["lower"]],
      upper = interval[["upper"]],
      posterior_width = interval[["width"]],
      width_success = interval[["width"]] <= config$width_target,
      covered = interval[["lower"]] <= true_estimand &&
        true_estimand <= interval[["upper"]],
      squared_error = (posterior_mean - true_estimand)^2,
      effective_sample_size = effective_sample_size,
      ess_fraction = effective_sample_size / config$particles,
      resampled = effective_sample_size <
        config$resampling_fraction * config$particles
    )

    if (effective_sample_size <
        config$resampling_fraction * config$particles) {
      indices <- systematic_resample(weights)
      particles <- particles[indices, , drop = FALSE]
      weights <- rep(1 / config$particles, config$particles)
    }
  }

  list(
    rounds = do.call(rbind, rows[!vapply(rows, is.null, logical(1))]),
    curves = if (length(curve_rows) > 0L) {
      do.call(rbind, curve_rows)
    } else {
      NULL
    }
  )
}

summarise_results <- function(round_results, config) {
  completed <- subset(round_results, feasible)
  policies <- unique(completed$policy)
  summaries <- lapply(policies, function(policy_name) {
    data <- subset(completed, policy == policy_name)
    totals <- aggregate(total_n ~ replicate, data, sum)
    data.frame(
      policy = policy_name,
      completed_replicates =
        length(unique(data$replicate[data$round == config$rounds])),
      decisions = nrow(data),
      mean_total_n = mean(data$total_n),
      median_total_n = median(data$total_n),
      mean_effort_per_replicate = mean(totals$total_n),
      width_success_rate = mean(data$width_success),
      interval_coverage = mean(data$covered),
      rmse = sqrt(mean(data$squared_error)),
      mean_predicted_assurance = mean(data$predicted_assurance),
      mean_assurance_lower_bound =
        mean(data$assurance_lower_bound),
      mean_validated_assurance = mean(data$validated_assurance),
      design_validation_difference =
        mean(data$predicted_assurance - data$validated_assurance),
      validation_calibration_error =
        mean(data$width_success) - mean(data$validated_assurance),
      mean_posterior_width = mean(data$posterior_width),
      mean_effective_sample_size = mean(data$effective_sample_size),
      resampling_rate = mean(data$resampled)
    )
  })
  do.call(rbind, summaries)
}

summarise_adaptive_periods <- function(round_results) {
  data <- subset(round_results, policy == "adaptive" & feasible)
  summaries <- aggregate(
    cbind(
      total_n,
      posterior_width,
      width_success,
      covered,
      predicted_assurance,
      validated_assurance
    ) ~ high_volatility,
    data,
    mean
  )
  names(summaries) <- c(
    "high_volatility",
    "mean_total_n",
    "mean_posterior_width",
    "width_success_rate",
    "interval_coverage",
    "mean_predicted_assurance",
    "mean_validated_assurance"
  )
  summaries$decisions <- as.integer(table(data$high_volatility)[
    as.character(summaries$high_volatility)
  ])
  summaries
}

run_replicate <- function(replicate_id, config) {
  set.seed(config$seed + replicate_id)
  truth <- simulate_truth(config)
  policy_results <- lapply(c("adaptive", "fixed_total"), function(policy) {
    run_policy(
      truth,
      policy,
      replicate_id,
      config
    )
  })
  message(
    "Completed replicate ",
    replicate_id,
    " of ",
    config$replicates
  )
  list(
    rounds = do.call(rbind, lapply(policy_results, `[[`, "rounds")),
    curves = policy_results[[1]]$curves
  )
}

make_plots <- function(round_results, assurance_curves, config) {
  example <- subset(round_results, replicate == 1 & feasible)
  adaptive_example <- subset(example, policy == "adaptive")

  trajectory <- ggplot(
    adaptive_example,
    aes(x = round, y = posterior_mean)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#8DB7C7",
      alpha = 0.45
    ) +
    geom_line(color = "#176B87", linewidth = 0.8) +
    geom_point(color = "#176B87", size = 1.8) +
    geom_line(
      aes(y = true_estimand),
      color = "#222222",
      linewidth = 0.7,
      linetype = 2
    ) +
    labs(
      x = NULL,
      y = "Population-weighted prevalence",
      title = "Adaptive filtering trajectory",
      subtitle = "Dashed line shows the simulated truth"
    )

  sample_size <- ggplot(
    subset(round_results, feasible),
    aes(x = round, y = total_n, color = policy)
  ) +
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 0.9
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.4
    ) +
    annotate(
      "rect",
      xmin = 5.5,
      xmax = 7.5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.10,
      fill = "#D95F02"
    ) +
    scale_color_manual(
      values = c(adaptive = "#176B87", fixed_total = "#8C3B2A"),
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    labs(
      x = NULL,
      y = "Total sample size",
      color = NULL,
      title = "Mean sampling effort",
      subtitle = "Rounds 6 and 7 have threefold innovation SD"
    )

  widths <- ggplot(
    subset(round_results, feasible),
    aes(x = posterior_width, fill = policy)
  ) +
    geom_histogram(
      position = "identity",
      alpha = 0.55,
      bins = 25
    ) +
    geom_vline(
      xintercept = config$width_target,
      linetype = 2,
      linewidth = 0.7
    ) +
    scale_fill_manual(
      values = c(adaptive = "#176B87", fixed_total = "#8C3B2A"),
      labels = c(adaptive = "Adaptive", fixed_total = "Fixed total")
    ) +
    labs(
      x = "Realized posterior interval width",
      y = "Round count",
      fill = NULL,
      title = "Posterior precision"
    )

  curve_round <- which.max(config$volatility_multiplier)
  curve_data <- subset(
    assurance_curves,
    replicate == 1 & round == curve_round
  )
  assurance_curve <- ggplot(
    curve_data,
    aes(x = total)
  ) +
    geom_line(
      aes(y = assurance, color = "Raw assurance"),
      linewidth = 0.75
    ) +
    geom_point(
      aes(y = assurance, color = "Raw assurance"),
      size = 1.5
    ) +
    geom_step(
      aes(y = sustained_assurance, color = "Sustained assurance"),
      linewidth = 0.75
    ) +
    geom_step(
      aes(
        y = sustained_lower_bound,
        color = "Sustained simultaneous lower bound"
      ),
      linewidth = 0.75,
      linetype = 3
    ) +
    geom_hline(
      yintercept = config$assurance_target,
      linetype = 2,
      linewidth = 0.6
    ) +
    geom_vline(
      data = subset(curve_data, selected),
      aes(xintercept = total),
      linewidth = 0.6,
      color = "#222222"
    ) +
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
      color = NULL,
      title = paste("Assurance curve in round", curve_round),
      subtitle = "Vertical line shows the selected total"
    )

  (trajectory / sample_size / assurance_curve / widths) +
    plot_layout(heights = c(1.05, 1, 1, 1)) &
    theme_minimal(base_size = 11) &
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

main <- function() {
  start_time <- Sys.time()
  options <- parse_options(commandArgs(trailingOnly = TRUE))
  config <- make_config(options)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  replicate_ids <- seq_len(config$replicates)
  if (config$workers > 1L && .Platform$OS.type != "windows") {
    replicate_results <- parallel::mclapply(
      replicate_ids,
      run_replicate,
      config = config,
      mc.cores = config$workers,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    replicate_results <- lapply(
      replicate_ids,
      run_replicate,
      config = config
    )
  }
  failed <- vapply(replicate_results, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("One or more simulation replicates failed")
  }

  round_results <- do.call(
    rbind,
    lapply(replicate_results, `[[`, "rounds")
  )
  assurance_curves <- do.call(
    rbind,
    lapply(replicate_results, `[[`, "curves")
  )
  summary <- summarise_results(round_results, config)
  period_summary <- summarise_adaptive_periods(round_results)
  example <- subset(round_results, replicate == 1)
  plot <- make_plots(round_results, assurance_curves, config)

  write.csv(
    round_results,
    file.path(config$output_dir, "simulation_rounds.csv"),
    row.names = FALSE
  )
  write.csv(
    summary,
    file.path(config$output_dir, "simulation_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    period_summary,
    file.path(config$output_dir, "adaptive_period_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    example,
    file.path(config$output_dir, "example_trajectory.csv"),
    row.names = FALSE
  )
  write.csv(
    assurance_curves,
    file.path(config$output_dir, "assurance_curves.csv"),
    row.names = FALSE
  )
  write.csv(
    subset(assurance_curves, replicate == 1),
    file.path(config$output_dir, "example_assurance_curves.csv"),
    row.names = FALSE
  )
  ggsave(
    file.path(config$output_dir, "simulation_overview.png"),
    plot,
    width = 8,
    height = 12,
    dpi = 180
  )
  ggsave(
    file.path(config$output_dir, "simulation_overview.pdf"),
    plot,
    width = 8,
    height = 12
  )

  config_table <- data.frame(
    quantity = c(
      "strata",
      "stratum weights",
      "rounds",
      "particles",
      "design Monte Carlo draws",
      "validation Monte Carlo draws",
      "selection family confidence",
      "candidate coupling",
      "credible level",
      "width target",
      "assurance target",
      "initial prevalence means",
      "initial logit standard deviations",
      "process innovation standard deviations",
      "process innovation correlation",
      "high-volatility rounds",
      "high-volatility innovation SD multiplier",
      "candidate totals",
      "fixed-total comparator",
      "replicates",
      "workers",
      "wall time seconds",
      "seed"
    ),
    value = c(
      paste(config$strata, collapse = ", "),
      paste(config$stratum_weights, collapse = ", "),
      config$rounds,
      config$particles,
      config$design_mc,
      config$validation_mc,
      config$selection_confidence,
      "common random numbers",
      config$credibility,
      config$width_target,
      config$assurance_target,
      paste(round(inv_logit(config$initial_mean), 4), collapse = ", "),
      paste(
        round(sqrt(diag(config$initial_covariance)), 4),
        collapse = ", "
      ),
      paste(
        round(sqrt(diag(config$process_covariance)), 4),
        collapse = ", "
      ),
      round(
        config$process_covariance[1, 2] /
          sqrt(
            config$process_covariance[1, 1] *
              config$process_covariance[2, 2]
          ),
        4
      ),
      paste(which(config$volatility_multiplier > 1), collapse = ", "),
      max(config$volatility_multiplier),
      paste(range(config$candidate_totals), collapse = " to "),
      config$fixed_total,
      config$replicates,
      config$workers,
      round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
      config$seed
    )
  )
  write.csv(
    config_table,
    file.path(config$output_dir, "simulation_config.csv"),
    row.names = FALSE
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(config$output_dir, "session_info.txt")
  )

  print(summary, row.names = FALSE)
}

main()
