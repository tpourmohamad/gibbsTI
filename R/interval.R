#' Combine marginal quantile posteriors into a two-sided interval
#'
#' Internal helper. Takes the alpha/2 quantile of the lower-limit posterior and
#' the 1 - alpha/2 quantile of the upper-limit posterior, where alpha = 1 - conf_level.
#' Kept in the package namespace (not exported) because the C++ calibration code
#' looks it up via `Environment::namespace_env("gibbsTI")`.
#'
#' @param q_lower Posterior draws for the lower quantile parameter.
#' @param q_upper Posterior draws for the upper quantile parameter.
#' @param conf_level Confidence level in (0, 1).
#' @return A list with numeric elements `lower` and `upper`.
#' @noRd
compute_two_sided_interval <- function(q_lower, q_upper, conf_level) {
  alpha <- 1 - conf_level
  # Take the alpha/2 and 1 - alpha/2 quantiles of the respective marginal posteriors.
  lower_bound <- stats::quantile(q_lower, alpha / 2)
  upper_bound <- stats::quantile(q_upper, 1 - alpha / 2)

  list(lower = as.numeric(lower_bound), upper = as.numeric(upper_bound))
}

# Compute the reported tolerance limit(s) from the posterior draws.
#
# One-sided: `samples` are posterior draws of the content-level quantile
# (targeted via `tau_use` in gibbsTI()). The reported limit is read off at the
# 1 - alpha = `confidence` level of those draws -- the upper confidence bound for
# an upper limit, the lower confidence bound (1 - confidence) for a lower limit.
.compute_interval <- function(samples, eta_val, side, type, content, confidence) {
  if (side == "one") {
    if (type == "upper") {
      return(unname(stats::quantile(samples, confidence)))
    } else {
      return(unname(stats::quantile(samples, 1 - confidence)))
    }
  } else {
    return(compute_two_sided_interval(samples[, 1], samples[, 2], confidence))
  }
}
