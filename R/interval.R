#' @export
compute_two_sided_interval <- function(q_lower, q_upper, conf_level) {
  alpha <- 1 - conf_level
  # Common approach: take the alpha/2 and 1-alpha/2
  # of the respective marginal posteriors
  lower_bound <- quantile(q_lower, alpha / 2)
  upper_bound <- quantile(q_upper, 1 - alpha / 2)

  list(lower = as.numeric(lower_bound), upper = as.numeric(upper_bound))
}

.compute_interval <- function(samples, eta_val, side, type, content, confidence) {
  if (side == "one") {
    if (type == "upper") {
      return(unname(quantile(samples, content)))
    } else {
      return(unname(quantile(samples, 1 - content)))
    }
  } else {
    return(compute_two_sided_interval(samples[,1], samples[,2], confidence))
  }
}
