#' Evaluate the operating characteristics of an interim design
#' @param design An object returned by [interim_design()].
#' @param drift Mean of the final Z statistic, from -10 to 10. The mean Z at
#'   information fraction `t` is `drift * sqrt(t)`. Use zero under the null.
#' @return A list with overall rejection probability, expected information
#'   fraction, and first-crossing and cumulative rejection probabilities by look.
#'   Trials that never cross consume the final information fraction.
#' @export
#' @examples
#' operating_characteristics(interim_design(c(0.5, 1)), drift = 2.8)
operating_characteristics <- function(design, drift = 0) {
  check_design(design, drift)
  probabilities <- crossing_probabilities(design$information, design$looks$boundary_z,
    drift = drift)$first_crossing
  summarize_crossing(design, probabilities)
}

#' Simulate interim analyses under the canonical joint-normal model
#' @inheritParams operating_characteristics
#' @param n_sim Integer number of trials, from 2 to 1000000.
#' @param seed Optional integer random seed. If supplied, the caller's
#'   `.Random.seed` is restored when the function returns.
#' @return The same fields as [operating_characteristics()], plus `n_sim`,
#'   `mc_se` (plug-in Monte Carlo SE), and `mc_ci` (95 percent exact binomial CI).
#' @export
#' @examples
#' simulate_interim(interim_design(c(0.5, 1)), n_sim = 10000, seed = 17)
simulate_interim <- function(design, drift = 0, n_sim = 100000L, seed = NULL) {
  check_design(design, drift)
  check_scalar(n_sim, "n_sim", 2, 1000000)
  if (n_sim != trunc(n_sim)) stop("`n_sim` must be an integer.", call. = FALSE)
  if (!is.null(seed)) {
    check_scalar(seed, "seed", 0, .Machine$integer.max)
    if (seed != trunc(seed)) stop("`seed` must be an integer.", call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", previous_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  information <- design$information
  first_cross <- integer(n_sim)
  increments <- diff(c(0, information))
  brownian <- numeric(n_sim)
  for (look in seq_along(information)) {
    brownian <- brownian + stats::rnorm(n_sim, sd = sqrt(increments[look]))
    z <- brownian / sqrt(information[look]) + drift * sqrt(information[look])
    first_cross[first_cross == 0L & z >= design$looks$boundary_z[look]] <- look
  }
  counts <- tabulate(first_cross, nbins = length(information))
  result <- summarize_crossing(design, counts / n_sim)
  result$n_sim <- n_sim
  result$mc_se <- sqrt(result$rejection_probability * (1 - result$rejection_probability) / n_sim)
  result$mc_ci <- as.numeric(stats::binom.test(sum(counts), n_sim)$conf.int)
  names(result$mc_ci) <- c("lower", "upper")
  result
}

check_design <- function(design, drift) {
  if (!inherits(design, "interim_design")) {
    stop("`design` must be created with interim_design().", call. = FALSE)
  }
  check_scalar(drift, "drift", -10, 10)
}

summarize_crossing <- function(design, probabilities) {
  information <- design$information
  list(rejection_probability = sum(probabilities),
    expected_information = 1 - sum((1 - information) * probabilities),
    looks = data.frame(look = seq_along(information), information = information,
      first_crossing = probabilities, cumulative_rejection = cumsum(probabilities)))
}
