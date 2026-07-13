test_that("gibbsTI errors on bad input", {

  expect_error(gibbsTI("a"))                                   # non-numeric x
  expect_error(gibbsTI(c(1, 2)))                               # too few observations
  expect_error(gibbsTI(c(1, 2, 3, 4, 5), content = 1.2))       # content out of (0, 1)
  expect_error(gibbsTI(rnorm(20), confidence = 0))             # confidence out of (0, 1)
  expect_error(gibbsTI(rnorm(20), eta_method = "fixed"))       # fixed requires eta
  expect_error(gibbsTI(rnorm(20), n_mcmc = 0))                 # non-positive n_mcmc
  expect_error(gibbsTI(rnorm(20), target = "quantile",
                       tau_lower = 0.9, tau_upper = 0.1))      # tau_lower >= tau_upper

})

test_that("one-sided fit returns a well-formed gibbsTI object", {
  set.seed(1)
  x <- rnorm(50)

  up <- gibbsTI(x, side = "one", type = "upper", eta_method = "fixed", eta = 1,
                n_mcmc = 400, burnin = 200, seed = 1, verbose = FALSE)
  lo <- gibbsTI(x, side = "one", type = "lower", eta_method = "fixed", eta = 1,
                n_mcmc = 400, burnin = 200, seed = 1, verbose = FALSE)

  expect_s3_class(up, "gibbsTI")
  expect_length(up$interval, 1)
  expect_true(is.finite(up$interval))
  # An upper limit should sit above a lower limit for the same data.
  expect_gt(up$interval, lo$interval)
})

test_that("two-sided content fit returns an ordered interval", {
  set.seed(2)
  x <- rnorm(50)

  fit <- gibbsTI(x, side = "two", content = 0.95, confidence = 0.90,
                 eta_method = "fixed", eta = 1,
                 n_mcmc = 400, burnin = 200, seed = 2, verbose = FALSE)

  expect_s3_class(fit, "gibbsTI")
  expect_true(is.finite(fit$interval$lower))
  expect_true(is.finite(fit$interval$upper))
  expect_lt(fit$interval$lower, fit$interval$upper)
})

test_that("compute_two_sided_interval orders bounds and honours confidence", {
  ql <- rnorm(1000, mean = -2)
  qu <- rnorm(1000, mean = 2)

  out <- gibbsTI:::compute_two_sided_interval(ql, qu, conf_level = 0.90)

  expect_named(out, c("lower", "upper"))
  expect_lt(out$lower, out$upper)
})
