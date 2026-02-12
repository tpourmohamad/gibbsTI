#' Internal wrapper for GPC calibration
.calibrate_eta_gpc <- function(x, side, target, content, confidence,
                               tau, tau_lower, tau_upper, verbose, eta0,
                               max_iter, B, tol, c, gamma, prop_sd, w, m,
                               n_samps_boot, burn_in_boot) {
  alpha <- 1 - confidence

  if (side == "one") {
    # Calls calibrate_eta_cpp in src/one_sided.cpp
    out <- calibrate_eta_cpp(
      y = x,
      tau = tau,
      alpha = alpha,
      eta0 = eta0,
      max_iter = max_iter,
      B = B,
      tol = tol,
      c = c,
      gamma = gamma,
      n_samps_boot = n_samps_boot,
      burn_in_boot = burn_in_boot,
      w = w,
      m = m,
      verbose = verbose
    )
  } else if (target == "content") {
    # Calls calibrate_eta_joint_content_cpp in src/two_sided_content.cpp
    out <- calibrate_eta_joint_content_cpp(
      y = x,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      alpha = alpha,
      eta0 = eta0,
      max_iter = max_iter,
      B = B,
      tol = tol,
      c = c,
      gamma = gamma,
      prop_sd = prop_sd,
      n_samps_boot = n_samps_boot,
      burn_in_boot = burn_in_boot,
      verbose = verbose
    )
  } else {
    # Calls calibrate_eta_joint_quantile_cpp in src/two_sided_quantile.cpp
    out <- calibrate_eta_joint_quantile_cpp(
      y = x,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      alpha = alpha,
      eta0 = eta0,
      max_iter = max_iter,
      B = B,
      tol = tol,
      c = c,
      gamma = gamma,
      prop_sd = prop_sd,
      n_samps_boot = n_samps_boot,
      burn_in_boot = burn_in_boot,
      verbose = verbose
    )
  }

  return(out$final_eta)
}
#' Internal wrapper for final MCMC sampling
.run_sampler <- function(x, side, target, tau, tau_lower, tau_upper,
                         eta, n_mcmc, burnin, thin, w, m, prop_sd, verbose) { # <--- Added prop_sd

  n_post_burnin <- n_mcmc * thin
  total_iterations <- burnin + n_post_burnin

  if (side == "one") {
    x0 <- stats::median(x)
    draws <- slice_sample_1d_cpp(
      x0 = x0, y = x, tau = tau, eta = eta,
      w = w, m = m, n_samps = total_iterations,
      lower = -Inf, upper = Inf
    )

    post_burnin_draws <- draws[(burnin + 1):total_iterations]
    final_samples <- post_burnin_draws[seq(1, length(post_burnin_draws), by = thin)]
    return(final_samples)

  } else {
    # 2D Samplers (Metropolis-Hastings)
    if (target == "content") {
      draws_matrix <- sample_joint_cpp(
        y = x, tau_lower = tau_lower, tau_upper = tau_upper,
        eta = eta,
        n_samps = n_post_burnin,
        burn_in = burnin,
        prop_sd = prop_sd
      )
    } else {
      draws_matrix <- sample_joint_quantile_cpp(
        y = x, tau_lower = tau_lower, tau_upper = tau_upper,
        eta = eta,
        n_samps = n_post_burnin,
        burn_in = burnin,
        prop_sd = prop_sd
      )
    }

    # Apply thinning to the matrix rows
    final_samples <- draws_matrix[seq(1, nrow(draws_matrix), by = thin), , drop = FALSE]
    return(final_samples)
  }
}
