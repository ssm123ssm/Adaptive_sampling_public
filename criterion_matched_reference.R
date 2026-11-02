posterior_width <- function(successes, total, tau) {
  shape1 <- successes + 0.5
  shape2 <- total - successes + 0.5
  qbeta(1 - tau / 2, shape1, shape2) - qbeta(tau / 2, shape1, shape2)
}

lower_passing_boundary <- function(total, width_target, tau) {
  peak <- as.integer(floor(total / 2))
  if (posterior_width(0L, total, tau) > width_target) {
    return(-1L)
  }
  if (posterior_width(peak, total, tau) <= width_target) {
    return(peak)
  }
  low <- 0L
  high <- peak
  while (high - low > 1L) {
    middle <- (low + high) %/% 2L
    if (posterior_width(middle, total, tau) <= width_target) {
      low <- middle
    } else {
      high <- middle
    }
  }
  low
}

fixed_sample_assurance <- function(total, prevalence, width_target, tau) {
  lower <- lower_passing_boundary(total, width_target, tau)
  upper <- total - lower
  passing <- if (lower < 0L) 0 else pbinom(lower, total, prevalence)
  if (upper <= total) {
    passing <- passing +
      pbinom(upper - 1L, total, prevalence, lower.tail = FALSE)
  }
  min(passing, 1)
}

criterion_matched_total <- function(prevalence, width_target, tau,
                                    assurance_target, upper_multiple = 3) {
  variance_matched <- (2 * qnorm(1 - tau / 2) *
    sqrt(prevalence * (1 - prevalence)) / width_target)^2
  lower <- max(1L, as.integer(floor(0.5 * variance_matched)))
  upper <- as.integer(ceiling(upper_multiple * variance_matched)) + 50L
  totals <- seq.int(lower, upper)
  assurance <- vapply(
    totals,
    fixed_sample_assurance,
    numeric(1),
    prevalence = prevalence,
    width_target = width_target,
    tau = tau
  )
  sustained <- rev(cummin(rev(assurance)))
  feasible <- which(sustained >= assurance_target)
  if (length(feasible) == 0L) {
    stop("The criterion-matched reference did not reach the assurance target")
  }
  totals[[feasible[[1]]]]
}

criterion_matched_curve <- function(prevalences, width_target, tau,
                                    assurance_target, spacing = 0.0025) {
  grid <- seq(
    max(spacing, floor(min(prevalences) / spacing) * spacing),
    ceiling(max(prevalences) / spacing) * spacing,
    by = spacing
  )
  tabulated <- vapply(
    grid,
    criterion_matched_total,
    numeric(1),
    width_target = width_target,
    tau = tau,
    assurance_target = assurance_target
  )
  approx(grid, tabulated, xout = prevalences, rule = 2)$y
}

variance_matched_total <- function(prevalence, width_target, tau) {
  (2 * qnorm(1 - tau / 2) *
    sqrt(prevalence * (1 - prevalence)) / width_target)^2
}

add_criterion_matched_reference <- function(results, config) {
  results$criterion_matched_total <- criterion_matched_curve(
    results$true_estimand,
    width_target = config$width_target,
    tau = 1 - config$credibility,
    assurance_target = config$assurance_target
  )
  results$criterion_matched_ratio <-
    results$total_n / results$criterion_matched_total
  results
}
