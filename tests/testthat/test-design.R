test_that("one look reduces to a fixed-horizon one-sided test", {
  for (method in c("obf", "pocock", "nominal")) {
    design <- interim_design(1, spending = method)
    expect_equal(design$looks$boundary_z, qnorm(0.975), tolerance = 1e-10)
    expect_equal(operating_characteristics(design, 2.8)$rejection_probability,
      pnorm(qnorm(0.975) - 2.8, lower.tail = FALSE), tolerance = 1e-10)
    expect_equal(operating_characteristics(design)$expected_information, 1)
  }
})

test_that("alpha spending targets hold across regular and irregular schedules", {
  for (info in list(c(0.5, 0.75, 1), c(0.2, 0.55, 0.91, 1), seq(0.1, 1, by = 0.1))) {
    obf <- interim_design(info)
    pocock <- interim_design(info, spending = "pocock")
    expect_equal(obf$looks$cumulative_alpha, 2 * pnorm(-qnorm(0.9875) / sqrt(info)), tolerance = 1e-10)
    expect_equal(pocock$looks$cumulative_alpha, 0.025 * log1p(expm1(1) * info), tolerance = 1e-10)
    expect_gt(sum(interim_design(info, spending = "nominal")$looks$alpha_increment), 0.025)
  }
})

test_that("crossing probabilities agree with an independent joint normal CDF", {
  skip_if_not_installed("mvtnorm")
  times <- c(0.2, 0.47, 0.72, 0.93, 1)
  corr <- sqrt(outer(times, times, pmin) / outer(times, times, pmax))
  for (method in c("obf", "pocock", "nominal")) {
    design <- interim_design(times, spending = method)
    for (drift in c(0, 1.2, 2.8)) {
      expected <- 1 - as.numeric(mvtnorm::pmvnorm(upper = design$looks$boundary_z,
        mean = drift * sqrt(times), corr = corr, algorithm = mvtnorm::Miwa(512)))
      expect_equal(operating_characteristics(design, drift)$rejection_probability, expected, tolerance = 1e-5)
    }
  }
})

test_that("seeded simulation is reproducible and preserves the caller's seed", {
  design <- interim_design(c(0.5, 0.75, 1))
  set.seed(29)
  previous <- .Random.seed
  result <- simulate_interim(design, n_sim = 200000L, seed = 905L)
  expect_identical(.Random.seed, previous)
  expect_identical(result, simulate_interim(design, n_sim = 200000L, seed = 905L))
  expect_lt(abs(result$rejection_probability - 0.025), 3 * result$mc_se)
  expect_lte(result$mc_ci[[1L]], result$rejection_probability)
  expect_gte(result$mc_ci[[2L]], result$rejection_probability)
  expect_equal(result$expected_information,
    sum(result$looks$first_crossing * design$information) + (1 - result$rejection_probability))
})

test_that("unsupported inputs get actionable errors", {
  expect_error(interim_design(c(0.05, 1)), "information")
  expect_error(interim_design(c(0.5, 0.5, 1)), "information")
  expect_error(interim_design(c(0.5, NA, 1)), "information")
  expect_error(interim_design(c(0.5, 0.8)), "information")
  expect_error(interim_design(1, alpha = 0), "alpha")
  expect_error(operating_characteristics(list()), "design")
  expect_error(simulate_interim(interim_design(1), n_sim = 5.5), "integer")
})
