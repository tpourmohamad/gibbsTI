#' Internal wrapper for GPC calibration
.calibrate_eta_gpc <- function(x, side, target, content, confidence,
                               tau, tau_lower, tau_upper, verbose, eta0) {
  alpha <- 1 - confidence

  if (side == "one") {
    # Calls calibrate_eta_cpp in src/one_sided.cpp
    out <- calibrate_eta_cpp(
      y = x,
      tau = tau,
      alpha = alpha,
      eta0 = eta0,      # Passed from R user input
      verbose = verbose
    )
  } else if (target == "content") {
    # Calls calibrate_eta_joint_content_cpp in src/two_sided_content.cpp
    out <- calibrate_eta_joint_content_cpp(
      y = x,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      alpha = alpha,
      eta0 = eta0,      # Passed from R user input
      verbose = verbose
    )
  } else {
    # Calls calibrate_eta_joint_quantile_cpp in src/two_sided_quantile.cpp
    out <- calibrate_eta_joint_quantile_cpp(
      y = x,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      alpha = alpha,
      eta0 = eta0,      # Passed from R user input
      verbose = verbose
    )
  }

  return(out$final_eta)
}
#' Internal wrapper for final MCMC sampling
.run_sampler <- function(x, side, target, tau, tau_lower, tau_upper, eta, n_mcmc, burnin, verbose) {

  # REMOVED: The message block is now handled by the parent function gibbsTI()

  if (side == "one") {
    # 1D Slice Sampler for One-Sided
    x0 <- median(x)
    draws <- slice_sample_1d_cpp(
      x0 = x0,
      y = x,
      tau = tau,
      eta = eta,
      w = 0.1,
      m = 100,
      n_samps = n_mcmc + burnin,
      lower = -Inf,
      upper = Inf
    )

    # Return only the post-burn-in vector
    return(draws[(burnin + 1):(n_mcmc + burnin)])

  } else {
    # 2D Samplers for Two-Sided
    if (target == "content") {
      draws_matrix <- sample_joint_cpp(
        y = x,
        tau_lower = tau_lower,
        tau_upper = tau_upper,
        eta = eta,
        n_samps = n_mcmc,
        burn_in = burnin
      )
    } else {
      draws_matrix <- sample_joint_quantile_cpp(
        y = x,
        tau_lower = tau_lower,
        tau_upper = tau_upper,
        eta = eta,
        n_samps = n_mcmc,
        burn_in = burnin
      )
    }

    return(draws_matrix)
  }
}
