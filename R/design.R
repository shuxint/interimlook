#' Design one-sided interim efficacy analyses
#'
#' Calibrate upper Z boundaries so that the cumulative null probability of
#' crossing by each look equals the specified alpha-spending function.
#' @param information Strictly increasing information fractions ending at 1.
#'   Supports one to ten looks, a first look at or after 0.10, and at least
#'   0.05 between successive looks.
#' @param alpha One-sided overall Type I error target, between 0.001 and 0.10.
#' @param spending One of `"obf"` (Lan-DeMets O'Brien-Fleming type),
#'   `"pocock"` (Lan-DeMets Pocock type), or `"nominal"` (apply `alpha`
#'   at every look as a benchmark for inflation).
#' @return An `interim_design` object with settings and a `looks` data frame.
#'   `boundary_z` is the efficacy cutoff. `nominal_p` is its one-sided
#'   marginal p-value. `alpha_increment` is the unconditional null probability
#'   of crossing for the first time at the look. `cumulative_alpha` sums those
#'   probabilities through each look.
#' @details The canonical model uses `Cov(Z_i, Z_j) = sqrt(t_i/t_j)` for
#'   `i <= j`. The O'Brien-Fleming type spends
#'   `2 * pnorm(qnorm(alpha/2, lower.tail = FALSE) / sqrt(t), lower.tail = FALSE)`.
#'   The Pocock type spends `alpha * log(1 + (exp(1) - 1) * t)`.
#'   Calibration assumes a single hypothesis and efficacy-only stopping.
#' @references Lan, K. K. G. and DeMets, D. L. (1983). Discrete sequential
#'   boundaries for clinical trials. Biometrika, 70(3), 659-663.
#'   \doi{10.1093/biomet/70.3.659}.
#' @export
#' @examples
#' interim_design(c(0.50, 0.75, 1))$looks
interim_design <- function(information, alpha = 0.025,
                           spending = c("obf", "pocock", "nominal")) {
  if (!is.numeric(information) || !is.null(dim(information)) ||
      !length(information) || length(information) > 10L ||
      any(!is.finite(information)) || information[1] < 0.10 ||
      tail_value(information) != 1 || any(diff(information) < 0.05 - 1e-12)) {
    stop("`information` must have 1 to 10 increasing fractions, starting at >= 0.10, ending at 1, with gaps >= 0.05.", call. = FALSE)
  }
  check_scalar(alpha, "alpha", 0.001, 0.10)
  spending <- match.arg(spending)
  target <- switch(spending,
    obf = 2 * stats::pnorm(stats::qnorm(alpha / 2, lower.tail = FALSE) /
      sqrt(information), lower.tail = FALSE),
    pocock = alpha * log1p(expm1(1) * information),
    nominal = NULL
  )
  boundary <- if (spending == "nominal") rep(stats::qnorm(alpha, lower.tail = FALSE), length(information))
  result <- crossing_probabilities(information, boundary, target = target)
  structure(list(
    alpha = alpha, spending = spending, information = information,
    looks = data.frame(look = seq_along(information), information = information,
      boundary_z = result$boundary, nominal_p = stats::pnorm(result$boundary, lower.tail = FALSE),
      alpha_increment = result$first_crossing, cumulative_alpha = cumsum(result$first_crossing))
  ), class = "interim_design")
}

#' Compare alpha spending and repeated nominal testing
#' @inheritParams interim_design
#' @return A long data frame with one row per method and analysis look.
#' @export
#' @examples
#' compare_interim(c(0.5, 1))
compare_interim <- function(information, alpha = 0.025) {
  methods <- c("obf", "pocock", "nominal")
  do.call(rbind, lapply(methods, function(method) {
    data.frame(spending = method, interim_design(information, alpha, method)$looks)
  }))
}

tail_value <- function(x) x[length(x)]

check_scalar <- function(x, name, lower, upper) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a finite number between %s and %s.", name, lower, upper), call. = FALSE)
  }
}

# Track a weighted density for paths that have stayed below all bounds. The
# transition uses centered Z coordinates, so drift only shifts each boundary.
crossing_probabilities <- function(information, boundary = NULL, drift = 0, target = NULL) {
  rule <- statmod::gauss.quad(96L)
  first_crossing <- numeric(length(information))
  increments <- if (!is.null(target)) diff(c(0, target))
  if (is.null(boundary)) boundary <- numeric(length(information))
  for (look in seq_along(information)) {
    time <- information[look]
    if (look == 1L) {
      crossing <- function(cutoff) stats::pnorm(cutoff, lower.tail = FALSE)
      density <- stats::dnorm
    } else {
      rho <- sqrt(information[look - 1L] / time)
      sd <- sqrt((time - information[look - 1L]) / time)
      means <- rho * nodes
      crossing <- function(cutoff) sum(mass * stats::pnorm((cutoff - means) / sd, lower.tail = FALSE))
      density <- function(x) drop(stats::dnorm(outer(x, means, "-"), sd = sd) %*% mass)
    }
    if (!is.null(target)) {
      boundary[look] <- stats::uniroot(function(z) crossing(z) - increments[look],
        c(-10, 18), tol = 1e-12)$root
    }
    cutoff <- boundary[look] - drift * sqrt(time)
    first_crossing[look] <- crossing(cutoff)
    if (look < length(information)) {
      # Normal mass below -12 or above +12 is below 2e-33 per look.
      upper <- max(-12, min(12, cutoff))
      width <- (upper + 12) / 2
      next_nodes <- -12 + width * (rule$nodes + 1)
      mass <- width * rule$weights * density(next_nodes)
      nodes <- next_nodes
    }
  }
  list(boundary = boundary, first_crossing = first_crossing)
}
