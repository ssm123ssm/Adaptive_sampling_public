suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

source("criterion_matched_reference.R")

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) plogis(x)

parse_options <- function(args) {
  defaults <- list(
    replicates = 100L,
    particles = 1500L,
    design_mc = 800L,
    validation_mc = 4000L,
    selection_confidence = 0.95,
    workers = 1L,
    rounds = 34L,
    variant_round = 12L,
    transition_rounds = 4L,
    initial_reproduction_number = 1.25,
    variant_multiplier = 1.55,
    gamma = 0.50,
    truth_beta_rw_sd = 0.04,
    beta_rw_sd = 0.10,
    width_target = 0.05,
    assurance_target = 0.90,
    credibility = 0.90,
    seed = 20260801L,
    output_dir = file.path(
      "output",
      "sir_variant_simulation"
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

make_config <- function(options) {
  contact_matrix <- matrix(
    c(1.20, 0.35,
      0.25, 0.90),
    nrow = 2,
    byrow = TRUE
  )
  dominant_eigenvalue <- max(Re(eigen(contact_matrix)$values))
  initial_beta <- options$initial_reproduction_number *
    options$gamma / dominant_eigenvalue

  config <- list(
    strata = c("A", "B"),
    stratum_weights = c(0.60, 0.40),
    contact_matrix = contact_matrix,
    dominant_eigenvalue = dominant_eigenvalue,
    initial_s = c(0.982, 0.989),
    initial_i = c(0.008, 0.005),
    initial_r = c(0.010, 0.006),
    initial_i_logit_sd = 0.25,
    initial_r_logit_sd = 0.15,
    initial_beta_log_sd = 0.10,
    initial_beta = initial_beta,
    initial_reproduction_number = options$initial_reproduction_number,
    variant_multiplier = options$variant_multiplier,
    variant_round = as.integer(options$variant_round),
    transition_rounds = as.integer(options$transition_rounds),
    gamma = options$gamma,
    rounds = as.integer(options$rounds),
    integration_substeps = 10L,
    truth_beta_rw_sd = options$truth_beta_rw_sd,
    beta_rw_sd = options$beta_rw_sd,
    beta_bounds = c(0.05, 3.00),
    particles = as.integer(options$particles),
    resampling_fraction = 0.50,
    credibility = options$credibility,
    width_target = options$width_target,
    assurance_target = options$assurance_target,
    candidate_totals = c(
      seq.int(20L, 200L, by = 20L),
      seq.int(240L, 600L, by = 40L),
      seq.int(700L, 1600L, by = 100L)
    ),
    stratum_minimum = c(10L, 10L),
    design_mc = as.integer(options$design_mc),
    validation_mc = as.integer(options$validation_mc),
    selection_confidence = options$selection_confidence,
    replicates = as.integer(options$replicates),
    workers = as.integer(options$workers),
    seed = as.integer(options$seed),
    output_dir = options$output_dir
  )

  if (config$variant_multiplier <= 1) {
    stop("The variant multiplier must exceed one")
  }
  if (min(config$candidate_totals) != sum(config$stratum_minimum)) {
    stop("The smallest candidate total must equal the summed stratum floor")
  }
  if (config$variant_round + config$transition_rounds >= config$rounds) {
    stop("The transition window must end before the final round")
  }
  config
}

round_period <- function(round_index, config) {
  ifelse(
    round_index < config$variant_round,
    "before",
    ifelse(
      round_index <= config$variant_round + config$transition_rounds,
      "transition",
      "after"
    )
  )
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
  sorted_values <- values[ordering]
  cumulative <- cumsum(weights[ordering]) / sum(weights)
  tail_probability <- (1 - credibility) / 2
  lower <- sorted_values[[which(cumulative >= tail_probability - 1e-12)[1]]]
  upper <- sorted_values[[which(cumulative >= 1 - tail_probability - 1e-12)[1]]]
  c(
    mean = sum(values * weights) / sum(weights),
    lower = lower,
    upper = upper,
    width = upper - lower
  )
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

sir_step <- function(s, i, r, beta, config) {
  step_size <- 1 / config$integration_substeps
  contact <- config$contact_matrix
  for (substep in seq_len(config$integration_substeps)) {
    pressure_a <- contact[1, 1] * i[, 1] + contact[1, 2] * i[, 2]
    pressure_b <- contact[2, 1] * i[, 1] + contact[2, 2] * i[, 2]
    force <- cbind(beta * pressure_a, beta * pressure_b)
    new_infections <- s * (1 - exp(-force * step_size))
    new_recoveries <- i * (1 - exp(-config$gamma * step_size))
    s <- s - new_infections
    i <- i + new_infections - new_recoveries
    r <- r + new_recoveries
  }
  list(s = s, i = i, r = r)
}

simulate_truth <- function(config) {
  s <- matrix(config$initial_s, nrow = 1)
  i <- matrix(config$initial_i, nrow = 1)
  r <- matrix(config$initial_r, nrow = 1)
  beta <- config$initial_beta
  rows <- vector("list", config$rounds)

  for (round_index in seq_len(config$rounds)) {
    beta <- beta * exp(rnorm(1, 0, config$truth_beta_rw_sd))
    if (round_index == config$variant_round) {
      beta <- beta * config$variant_multiplier
    }
    beta <- min(max(beta, config$beta_bounds[[1]]), config$beta_bounds[[2]])
    state <- sir_step(s, i, r, beta, config)
    s <- state$s
    i <- state$i
    r <- state$r
    rows[[round_index]] <- data.frame(
      round = round_index,
      true_beta = beta,
      true_reproduction_number =
        beta * config$dominant_eigenvalue / config$gamma,
      true_i = sum(config$stratum_weights * i[1, ]),
      true_r_A = r[1, 1],
      true_r_B = r[1, 2],
      true_estimand = sum(config$stratum_weights * r[1, ])
    )
  }
  do.call(rbind, rows)
}

initialise_particles <- function(config) {
  particle_count <- config$particles
  i <- sapply(config$initial_i, function(value) {
    inv_logit(logit(value) + rnorm(particle_count, 0, config$initial_i_logit_sd))
  })
  r <- sapply(config$initial_r, function(value) {
    inv_logit(logit(value) + rnorm(particle_count, 0, config$initial_r_logit_sd))
  })
  s <- 1 - i - r
  if (any(s <= 0)) {
    stop("Initial particle states left no susceptible population")
  }
  beta <- exp(
    log(config$initial_beta) +
      rnorm(particle_count, 0, config$initial_beta_log_sd)
  )
  list(s = s, i = i, r = r, beta = beta)
}

propagate_particles <- function(particles, config) {
  particles$beta <- exp(
    log(particles$beta) +
      rnorm(length(particles$beta), 0, config$beta_rw_sd)
  )
  particles$beta <- pmin(
    pmax(particles$beta, config$beta_bounds[[1]]),
    config$beta_bounds[[2]]
  )
  state <- sir_step(
    particles$s,
    particles$i,
    particles$r,
    particles$beta,
    config
  )
  particles$s <- state$s
  particles$i <- state$i
  particles$r <- state$r
  particles
}

build_greedy_family <- function(probabilities, predictive_weights, config) {
  coefficients <- config$stratum_weights^2 *
    colSums(probabilities * (1 - probabilities) * predictive_weights)
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

make_predictive_mc_stream <- function(
  probabilities,
  predictive_weights,
  maximum_allocation,
  draws
) {
  ancestors <- sample.int(
    nrow(probabilities),
    draws,
    replace = TRUE,
    prob = predictive_weights
  )
  cumulative_counts <- lapply(
    seq_along(maximum_allocation),
    function(stratum) {
      maximum <- maximum_allocation[[stratum]]
      uniforms <- matrix(runif(draws * maximum), nrow = draws)
      successes <- uniforms <= probabilities[ancestors, stratum]
      t(apply(successes, 1, cumsum))
    }
  )
  list(cumulative_counts = cumulative_counts, draws = draws)
}

counts_from_mc_stream <- function(stream, allocation) {
  counts <- matrix(0L, nrow = stream$draws, ncol = length(allocation))
  for (stratum in seq_along(allocation)) {
    counts[, stratum] <-
      stream$cumulative_counts[[stratum]][, allocation[[stratum]]]
  }
  counts
}

hypothetical_widths <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation,
  prospective_counts,
  credibility
) {
  clipped <- pmin(pmax(probabilities, 1e-10), 1 - 1e-10)
  log_probabilities <- log(clipped)
  log_complements <- log1p(-clipped)
  log_weight_matrix <- matrix(
    log(pmax(predictive_weights, 1e-300)),
    nrow = nrow(probabilities),
    ncol = nrow(prospective_counts)
  )
  for (stratum in seq_along(allocation)) {
    log_weight_matrix <- log_weight_matrix +
      tcrossprod(
        log_probabilities[, stratum],
        prospective_counts[, stratum]
      ) +
      tcrossprod(
        log_complements[, stratum],
        allocation[[stratum]] - prospective_counts[, stratum]
      )
  }
  maxima <- apply(log_weight_matrix, 2, max)
  posterior_weights <- exp(sweep(log_weight_matrix, 2, maxima, FUN = "-"))
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

estimate_assurance <- function(
  probabilities,
  estimands,
  predictive_weights,
  allocation,
  config,
  mc_draws,
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
      allocation <- allocation_family[index, ]
      estimate_assurance(
        probabilities,
        estimands,
        predictive_weights,
        allocation,
        config,
        mc_draws = config$design_mc,
        prospective_counts = counts_from_mc_stream(design_stream, allocation),
        comparisons = candidate_count
      )
    },
    numeric(4)
  ))

  raw_assurance <- evaluations[, "assurance"]
  lower_bounds <- evaluations[, "assurance_lower_bound"]
  sustained_assurance <- rev(cummin(rev(raw_assurance)))
  sustained_lower_bound <- rev(cummin(rev(lower_bounds)))
  feasible <- which(sustained_lower_bound >= config$assurance_target)
  selected <- if (length(feasible) == 0L) candidate_count else feasible[[1]]

  list(
    feasible = length(feasible) > 0L,
    index = selected,
    total = config$candidate_totals[[selected]],
    allocation = allocation_family[selected, ],
    design_assurance = raw_assurance[[selected]],
    sustained_assurance = sustained_assurance[[selected]],
    sustained_lower_bound = sustained_lower_bound[[selected]],
    nonmonotone = any(diff(raw_assurance) < 0),
    raw_sustained_gap = max(raw_assurance - sustained_assurance)
  )
}

run_replicate <- function(replicate_id, config) {
  set.seed(config$seed + replicate_id * 7919L)
  truth <- simulate_truth(config)
  particles <- initialise_particles(config)
  weights <- rep(1 / config$particles, config$particles)
  z_value <- qnorm(1 - (1 - config$credibility) / 2)
  rows <- vector("list", config$rounds)

  for (round_index in seq_len(config$rounds)) {
    particles <- propagate_particles(particles, config)
    estimands <- as.vector(particles$r %*% config$stratum_weights)
    predictive_estimand <- weighted_interval(
      estimands,
      weights,
      config$credibility
    )
    predictive_beta <- weighted_interval(
      particles$beta,
      weights,
      config$credibility
    )

    allocation_family <- build_greedy_family(particles$r, weights, config)
    design <- select_adaptive_design(
      particles$r,
      estimands,
      weights,
      allocation_family,
      config
    )
    validation <- estimate_assurance(
      particles$r,
      estimands,
      weights,
      design$allocation,
      config,
      mc_draws = config$validation_mc
    )

    truth_row <- truth[round_index, ]
    true_probabilities <- c(truth_row$true_r_A, truth_row$true_r_B)
    counts <- vapply(
      seq_along(design$allocation),
      function(stratum) {
        rbinom(
          1,
          design$allocation[[stratum]],
          true_probabilities[[stratum]]
        )
      },
      integer(1)
    )

    weights <- normalise_log_weights(
      log(pmax(weights, 1e-300)) +
        particle_log_likelihood(particles$r, counts, design$allocation)
    )
    filtering_estimand <- weighted_interval(
      estimands,
      weights,
      config$credibility
    )
    filtering_beta <- weighted_interval(
      particles$beta,
      weights,
      config$credibility
    )
    effective_sample_size <- 1 / sum(weights^2)
    oracle_fixed_total <- (2 * z_value *
      sqrt(truth_row$true_estimand * (1 - truth_row$true_estimand)) /
      config$width_target)^2

    rows[[round_index]] <- data.frame(
      replicate = replicate_id,
      round = round_index,
      period = round_period(round_index, config),
      true_beta = truth_row$true_beta,
      true_reproduction_number = truth_row$true_reproduction_number,
      true_i = truth_row$true_i,
      true_estimand = truth_row$true_estimand,
      predictive_mean = predictive_estimand[["mean"]],
      predictive_lower = predictive_estimand[["lower"]],
      predictive_upper = predictive_estimand[["upper"]],
      predictive_width = predictive_estimand[["width"]],
      predictive_covered =
        predictive_estimand[["lower"]] <= truth_row$true_estimand &&
        truth_row$true_estimand <= predictive_estimand[["upper"]],
      filtering_mean = filtering_estimand[["mean"]],
      filtering_lower = filtering_estimand[["lower"]],
      filtering_upper = filtering_estimand[["upper"]],
      filtering_width = filtering_estimand[["width"]],
      filtering_covered =
        filtering_estimand[["lower"]] <= truth_row$true_estimand &&
        truth_row$true_estimand <= filtering_estimand[["upper"]],
      width_success = filtering_estimand[["width"]] <= config$width_target,
      beta_predictive_mean = predictive_beta[["mean"]],
      beta_filtering_mean = filtering_beta[["mean"]],
      beta_filtering_lower = filtering_beta[["lower"]],
      beta_filtering_upper = filtering_beta[["upper"]],
      beta_filtering_width = filtering_beta[["width"]],
      beta_filtering_covered =
        filtering_beta[["lower"]] <= truth_row$true_beta &&
        truth_row$true_beta <= filtering_beta[["upper"]],
      beta_relative_error =
        abs(filtering_beta[["mean"]] - truth_row$true_beta) /
          truth_row$true_beta,
      total_n = design$total,
      n_A = design$allocation[[1]],
      n_B = design$allocation[[2]],
      y_A = counts[[1]],
      y_B = counts[[2]],
      oracle_fixed_total = oracle_fixed_total,
      effort_ratio = design$total / oracle_fixed_total,
      feasible = design$feasible,
      infeasible = !design$feasible,
      design_assurance = design$design_assurance,
      sustained_lower_bound = design$sustained_lower_bound,
      validation_assurance = validation[["assurance"]],
      validation_mc_se = validation[["mc_se"]],
      design_validation_gap =
        abs(design$design_assurance - validation[["assurance"]]),
      nonmonotone = design$nonmonotone,
      raw_sustained_gap = design$raw_sustained_gap,
      effective_sample_size = effective_sample_size,
      ess_fraction = effective_sample_size / config$particles,
      resampled = effective_sample_size <
        config$resampling_fraction * config$particles
    )

    if (effective_sample_size <
        config$resampling_fraction * config$particles) {
      indices <- systematic_resample(weights)
      particles$s <- particles$s[indices, , drop = FALSE]
      particles$i <- particles$i[indices, , drop = FALSE]
      particles$r <- particles$r[indices, , drop = FALSE]
      particles$beta <- particles$beta[indices]
      weights <- rep(1 / config$particles, config$particles)
    }
  }
  do.call(rbind, rows)
}

first_sustained_round <- function(rounds, condition, start_round, run = 2L) {
  eligible <- which(rounds >= start_round)
  for (index in eligible) {
    end <- index + run - 1L
    if (end <= length(condition) && all(condition[index:end])) {
      return(rounds[[index]])
    }
  }
  NA_integer_
}

make_recovery_summary <- function(results, config) {
  by_replicate <- split(results, results$replicate)
  recovery <- do.call(rbind, lapply(by_replicate, function(data) {
    data <- data[order(data$round), ]
    beta_close <- data$beta_relative_error <= 0.15
    trajectory_close <-
      abs(data$filtering_mean - data$true_estimand) <= 0.015 &
        data$filtering_covered
    departures <- which(
      data$round >= config$variant_round &
        (abs(data$predictive_mean - data$true_estimand) > 0.010 |
           !data$predictive_covered)
    )
    departure_round <- if (length(departures) == 0L) {
      NA_integer_
    } else {
      data$round[departures[[1]]]
    }
    trajectory_recovery_round <- if (is.na(departure_round)) {
      NA_integer_
    } else {
      first_sustained_round(
        data$round,
        trajectory_close,
        departure_round + 1L
      )
    }
    data.frame(
      replicate = data$replicate[[1]],
      beta_recovery_round = first_sustained_round(
        data$round,
        beta_close,
        config$variant_round
      ),
      trajectory_departure_round = departure_round,
      trajectory_recovery_round = trajectory_recovery_round
    )
  }))
  recovery$beta_recovery_lag <-
    recovery$beta_recovery_round - config$variant_round
  recovery$trajectory_recovery_lag <-
    recovery$trajectory_recovery_round - recovery$trajectory_departure_round
  recovery
}

replicate_measure <- function(results, column, subset = NULL) {
  data <- if (is.null(subset)) results else results[subset, , drop = FALSE]
  per_replicate <- tapply(data[[column]], data$replicate, mean)
  as.numeric(per_replicate)
}

monte_carlo_interval <- function(values, bounded = c("none", "unit")) {
  bounded <- match.arg(bounded)
  estimate <- mean(values, na.rm = TRUE)
  usable <- values[is.finite(values)]
  standard_error <- stats::sd(usable) / sqrt(length(usable))
  lower <- estimate - 1.96 * standard_error
  upper <- estimate + 1.96 * standard_error
  if (bounded == "unit") {
    lower <- max(lower, 0)
    upper <- min(upper, 1)
  }
  c(estimate = estimate, lower = lower, upper = upper)
}

make_component_summary <- function(results, recovery, config) {
  periods <- c("before", "transition", "after")
  labels <- c(
    before = "Before the jump",
    transition = "Transition window",
    after = "After the window"
  )
  rows <- list()
  add <- function(component, measure, values, bounded = "none") {
    interval <- monte_carlo_interval(values, bounded)
    rows[[length(rows) + 1L]] <<- data.frame(
      component = component,
      measure = measure,
      estimate = interval[["estimate"]],
      lower = interval[["lower"]],
      upper = interval[["upper"]]
    )
  }

  for (period in periods) {
    subset <- results$period == period
    add(
      "Interval coverage",
      paste0("Predictive coverage, ", labels[[period]]),
      replicate_measure(results, "predictive_covered", subset),
      "unit"
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Interval coverage",
      paste0("Posterior coverage, ", labels[[period]]),
      replicate_measure(results, "filtering_covered", subset),
      "unit"
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Transmission recovery",
      paste0("Transmission-rate coverage, ", labels[[period]]),
      replicate_measure(results, "beta_filtering_covered", subset),
      "unit"
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Transmission recovery",
      paste0("Transmission-rate relative error, ", labels[[period]]),
      replicate_measure(results, "beta_relative_error", subset)
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Design response",
      paste0("Mean selected total, ", labels[[period]]),
      replicate_measure(results, "total_n", subset)
    )
  }
  for (period in periods) {
    subset <- results$period == period &
      results$total_n > min(config$candidate_totals)
    add(
      "Design response",
      paste0("Effort relative to fixed-sample rule, ", labels[[period]]),
      replicate_measure(results, "criterion_matched_ratio", subset)
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Assurance",
      paste0("Realized width success, ", labels[[period]]),
      replicate_measure(results, "width_success", subset),
      "unit"
    )
  }
  for (period in periods) {
    subset <- results$period == period
    add(
      "Assurance",
      paste0("Validated assurance, ", labels[[period]]),
      replicate_measure(results, "validation_assurance", subset),
      "unit"
    )
  }
  add(
    "Recovery timing",
    "Transmission-rate recovery lag in rounds",
    recovery$beta_recovery_lag[is.finite(recovery$beta_recovery_lag)]
  )
  add(
    "Recovery timing",
    "Trajectory recovery lag in rounds",
    recovery$trajectory_recovery_lag[
      is.finite(recovery$trajectory_recovery_lag)
    ]
  )
  add(
    "Numerical",
    "Mean effective-sample-size fraction",
    replicate_measure(results, "ess_fraction"),
    "unit"
  )
  add(
    "Numerical",
    "Design-validation assurance error",
    replicate_measure(results, "design_validation_gap")
  )
  add(
    "Numerical",
    "Fraction of nonmonotone round curves",
    replicate_measure(results, "nonmonotone"),
    "unit"
  )
  add(
    "Numerical",
    "Fraction of infeasible decisions",
    replicate_measure(results, "infeasible"),
    "unit"
  )
  do.call(rbind, rows)
}

write_component_table <- function(summary_frame, config, path) {
  formatted <- summary_frame
  numeric_columns <- c("estimate", "lower", "upper")
  digits <- ifelse(
    grepl("Mean selected total|lag in rounds", formatted$measure),
    1,
    3
  )
  for (column in numeric_columns) {
    formatted[[column]] <- mapply(
      function(value, d) formatC(value, format = "f", digits = d),
      formatted[[column]],
      digits
    )
  }
  body <- paste0(
    formatted$component, " & ", formatted$measure, " & ",
    formatted$estimate, " & (", formatted$lower, ", ", formatted$upper, ") \\\\"
  )
  duplicated_component <- duplicated(formatted$component)
  body[duplicated_component] <- sub("^[^&]*&", " &", body[duplicated_component])

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Stratified SIR simulation with an unanticipated change in transmission. ",
      "Summaries are across ", config$replicates,
      " replicated trajectories. Intervals are normal-approximation 95\\% ",
      "Monte Carlo intervals for the replicate-level mean.}"
    ),
    "\\label{tab:sir_variant_summary}",
    "\\footnotesize",
    "\\begin{tabular}{llrr}",
    "\\toprule",
    "Component & Measure & Estimate & 95\\% interval \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, path)
}

round_summary_frame <- function(results) {
  columns <- c(
    "true_estimand", "true_beta", "true_reproduction_number",
    "predictive_mean", "predictive_lower", "predictive_upper",
    "predictive_width", "predictive_covered",
    "filtering_mean", "filtering_lower", "filtering_upper",
    "filtering_width", "filtering_covered", "width_success",
    "beta_filtering_mean", "beta_filtering_lower", "beta_filtering_upper",
    "beta_filtering_covered", "beta_relative_error",
    "total_n", "oracle_fixed_total", "effort_ratio",
    "criterion_matched_total", "criterion_matched_ratio",
    "design_assurance", "validation_assurance", "ess_fraction"
  )
  aggregated <- aggregate(
    results[, columns],
    by = list(round = results$round),
    FUN = mean
  )
  aggregated
}

separate_panel_rows <- function(figure) {
  figure & theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(6, 8, 16, 8)
  )
}

make_dynamics_plot <- function(results, config) {
  summary_frame <- round_summary_frame(results)
  band <- data.frame(
    xmin = config$variant_round - 0.5,
    xmax = config$variant_round + config$transition_rounds + 0.5
  )
  shade <- function(plot, ymin = -Inf, ymax = Inf) {
    plot + geom_rect(
      data = band,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#f4a261",
      alpha = 0.15
    )
  }

  p_trajectory <- shade(ggplot(summary_frame, aes(x = round))) +
    geom_ribbon(
      aes(ymin = predictive_lower, ymax = predictive_upper),
      fill = "grey65",
      alpha = 0.55
    ) +
    geom_ribbon(
      aes(ymin = filtering_lower, ymax = filtering_upper),
      fill = "#4c78a8",
      alpha = 0.35
    ) +
    geom_line(aes(y = true_estimand), colour = "black", linewidth = 0.8) +
    geom_line(aes(y = filtering_mean), colour = "#2c6aa0", linewidth = 0.8) +
    labs(
      x = "Sampling round",
      y = "Seroprevalence",
      title = "Estimand recovery"
    ) +
    theme_minimal(base_size = 10)

  p_beta <- shade(ggplot(summary_frame, aes(x = round))) +
    geom_ribbon(
      aes(ymin = beta_filtering_lower, ymax = beta_filtering_upper),
      fill = "#59a14f",
      alpha = 0.30
    ) +
    geom_line(aes(y = true_beta), colour = "black", linewidth = 0.8) +
    geom_line(
      aes(y = beta_filtering_mean),
      colour = "#3b7f36",
      linewidth = 0.8
    ) +
    labs(
      x = "Sampling round",
      y = expression(beta[t]),
      title = "Transmission-rate recovery"
    ) +
    theme_minimal(base_size = 10)

  width_frame <- rbind(
    data.frame(
      round = summary_frame$round,
      series = "One-step predictive",
      width = summary_frame$predictive_width
    ),
    data.frame(
      round = summary_frame$round,
      series = "Posterior",
      width = summary_frame$filtering_width
    )
  )
  p_width <- shade(ggplot(width_frame, aes(x = round, y = width))) +
    geom_hline(
      yintercept = config$width_target,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_line(aes(colour = series), linewidth = 0.8) +
    scale_colour_manual(
      values = c("One-step predictive" = "grey35", "Posterior" = "#2c6aa0")
    ) +
    labs(
      x = "Sampling round",
      y = "Interval width",
      colour = NULL,
      title = "Interval width around the change"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  effort_frame <- rbind(
    data.frame(
      round = summary_frame$round,
      series = "Adaptive policy",
      total = summary_frame$total_n
    ),
    data.frame(
      round = summary_frame$round,
      series = "Fixed-sample reference",
      total = summary_frame$criterion_matched_total
    )
  )
  p_effort <- shade(
    ggplot(effort_frame, aes(x = round, y = total)),
    ymin = 0.8 * min(effort_frame$total),
    ymax = 1.2 * max(effort_frame$total)
  ) +
    geom_line(aes(colour = series), linewidth = 0.8) +
    scale_colour_manual(
      values = c(
        "Adaptive policy" = "#b24a33",
        "Fixed-sample reference" = "grey45"
      )
    ) +
    scale_y_continuous(trans = "log10") +
    labs(
      x = "Sampling round",
      y = "Selected total",
      colour = NULL,
      title = "Sampling effort"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  separate_panel_rows((p_trajectory | p_beta) / (p_width | p_effort))
}

make_calibration_plot <- function(results, config) {
  summary_frame <- round_summary_frame(results)
  band <- data.frame(
    xmin = config$variant_round - 0.5,
    xmax = config$variant_round + config$transition_rounds + 0.5
  )
  shade <- function(plot) {
    plot + geom_rect(
      data = band,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#f4a261",
      alpha = 0.15
    )
  }

  coverage_frame <- rbind(
    data.frame(
      round = summary_frame$round,
      series = "One-step predictive",
      value = summary_frame$predictive_covered
    ),
    data.frame(
      round = summary_frame$round,
      series = "Posterior",
      value = summary_frame$filtering_covered
    )
  )
  p_coverage <- shade(ggplot(coverage_frame, aes(x = round, y = value))) +
    geom_hline(
      yintercept = config$credibility,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_line(aes(colour = series), linewidth = 0.8) +
    scale_colour_manual(
      values = c("One-step predictive" = "grey35", "Posterior" = "#2c6aa0")
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Sampling round",
      y = "Coverage",
      colour = NULL,
      title = "Interval coverage"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  assurance_frame <- rbind(
    data.frame(
      round = summary_frame$round,
      series = "Validated assurance",
      value = summary_frame$validation_assurance
    ),
    data.frame(
      round = summary_frame$round,
      series = "Realized width success",
      value = summary_frame$width_success
    )
  )
  p_assurance <- shade(ggplot(assurance_frame, aes(x = round, y = value))) +
    geom_hline(
      yintercept = config$assurance_target,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_line(aes(colour = series), linewidth = 0.8) +
    scale_colour_manual(
      values = c(
        "Validated assurance" = "#7b5aa6",
        "Realized width success" = "#b24a33"
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Sampling round",
      y = "Probability",
      colour = NULL,
      title = "Assurance against realized width success"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")

  p_beta_error <- shade(
    ggplot(summary_frame, aes(x = round, y = beta_relative_error))
  ) +
    geom_line(colour = "#3b7f36", linewidth = 0.8) +
    labs(
      x = "Sampling round",
      y = "Relative error",
      title = "Transmission-rate relative error"
    ) +
    theme_minimal(base_size = 10)

  p_ess <- shade(ggplot(summary_frame, aes(x = round, y = ess_fraction))) +
    geom_hline(
      yintercept = config$resampling_fraction,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_line(colour = "#8a6d3b", linewidth = 0.8) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Sampling round",
      y = "ESS fraction",
      title = "Effective sample size"
    ) +
    theme_minimal(base_size = 10)

  separate_panel_rows((p_coverage | p_assurance) / (p_beta_error | p_ess))
}

main <- function() {
  options <- parse_options(commandArgs(trailingOnly = TRUE))
  config <- make_config(options)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  started <- Sys.time()
  replicate_ids <- seq_len(config$replicates)
  results_list <- if (config$workers > 1L) {
    mclapply(replicate_ids, run_replicate, config = config,
             mc.cores = config$workers)
  } else {
    lapply(replicate_ids, run_replicate, config = config)
  }
  failed <- vapply(results_list, function(x) !is.data.frame(x), logical(1))
  if (any(failed)) {
    stop("Replicates failed: ", paste(which(failed), collapse = ", "))
  }
  results <- do.call(rbind, results_list)
  results <- add_criterion_matched_reference(results, config)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  recovery <- make_recovery_summary(results, config)
  component_summary <- make_component_summary(results, recovery, config)
  rounds_frame <- round_summary_frame(results)

  write.csv(
    results,
    file.path(config$output_dir, "sir_variant_round_results.csv"),
    row.names = FALSE
  )
  write.csv(
    rounds_frame,
    file.path(config$output_dir, "sir_variant_round_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    recovery,
    file.path(config$output_dir, "sir_variant_recovery_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    component_summary,
    file.path(config$output_dir, "sir_variant_component_summary.csv"),
    row.names = FALSE
  )
  write_component_table(
    component_summary,
    config,
    file.path(config$output_dir, "sir_variant_component_summary_table.tex")
  )

  dynamics_plot <- make_dynamics_plot(results, config)
  calibration_plot <- make_calibration_plot(results, config)
  for (extension in c("pdf", "png")) {
    ggsave(
      file.path(config$output_dir, paste0("sir_variant_dynamics.", extension)),
      dynamics_plot,
      width = 10,
      height = 7,
      dpi = 200
    )
    ggsave(
      file.path(
        config$output_dir,
        paste0("sir_variant_calibration.", extension)
      ),
      calibration_plot,
      width = 10,
      height = 7,
      dpi = 200
    )
  }

  configuration <- data.frame(
    parameter = c(
      "replicates", "particles", "design draws", "validation draws",
      "selection confidence", "rounds", "variant round", "transition rounds",
      "initial reproduction number", "variant multiplier", "gamma",
      "truth beta random-walk SD", "filter beta random-walk SD",
      "credibility", "width target", "assurance target",
      "candidate totals", "stratum floor", "seed", "runtime minutes"
    ),
    value = c(
      config$replicates, config$particles, config$design_mc,
      config$validation_mc, config$selection_confidence, config$rounds,
      config$variant_round, config$transition_rounds,
      config$initial_reproduction_number, config$variant_multiplier,
      config$gamma, config$truth_beta_rw_sd, config$beta_rw_sd,
      config$credibility, config$width_target, config$assurance_target,
      paste(range(config$candidate_totals), collapse = " to "),
      paste(config$stratum_minimum, collapse = ", "),
      config$seed, round(elapsed, 2)
    )
  )
  write.csv(
    configuration,
    file.path(config$output_dir, "sir_variant_config.csv"),
    row.names = FALSE
  )

  cat("Completed", config$replicates, "replicates in",
      sprintf("%.1f", elapsed), "minutes\n\n")
  print(component_summary, row.names = FALSE, digits = 3)
}

main()
