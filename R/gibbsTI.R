#' Fit a Gibbs Posterior Tolerance Interval
#'
#' @param x Numeric vector of observations.
#' @param side Direction of the interval: "one" (one-sided) or "two" (two-sided centered).
#' @param type Type of one-sided interval: "upper" (bound above) or "lower" (bound below).
#' @param target The calibration target: "content" (default for TIs) or "quantile".
#' @param content Population content level (e.g., 0.95).
#' @param confidence Confidence level (e.g., 0.90).
#' @param tau Target quantile (required for one-sided quantile target).
#' @param tau_lower Lower quantile (for two-sided quantile target).
#' @param tau_upper Upper quantile (for two-sided quantile target).
#' @param calibrate Logical; whether to use GPC calibration.
#' @param eta_method "GPC" (automated) or "fixed".
#' @param eta Fixed value of eta if calibrate = FALSE.
#' @param n_mcmc Number of iterations for final sampling.
#' @param burnin Number of burn-in iterations.
#' @param seed Random seed.
#' @param verbose Logical; if TRUE, prints progress messages and iterative calibration updates to the console.
#'
#' @return An object of class \code{"gibbsTI"}.
#' @export
gibbsTI <- function(
    x,
    side = c("one", "two"),
    type = c("upper", "lower"),
    target = c("content", "quantile"),
    content = 0.95,
    confidence = 0.90,
    tau = NULL,
    tau_lower = NULL,
    tau_upper = NULL,
    calibrate = TRUE,
    eta_method = c("GPC", "fixed"),
    eta = NULL,
    n_mcmc = 5000,
    burnin = 1000,
    seed = NULL,
    verbose = TRUE
) {

  # 1. Match Arguments
  side <- match.arg(side)
  type <- match.arg(type)
  target <- match.arg(target)
  eta_method <- match.arg(eta_method)

  # 2. Logic to set tau values
  if (target == "content") {
    if (side == "one") {
      tau_use <- if(type == "upper") content else (1 - content)
    } else {
      alpha_content <- 1 - content
      tau_lower <- alpha_content / 2
      tau_upper <- 1 - (alpha_content / 2)
      tau_use <- NULL
    }
  } else {
    # If target is 'quantile', prioritize user-supplied taus,
    # but fallback to content-based tails if user left them NULL
    tau_use <- if(!is.null(tau)) tau else (if(type == "upper") content else (1-content))

    if (side == "two") {
      if(is.null(tau_lower)) tau_lower <- (1 - content) / 2
      if(is.null(tau_upper)) tau_upper <- 1 - (1 - content) / 2
    }
  }

  # 3. Validation
  .validate_inputs(x, content, confidence, target, tau_use, tau_lower, tau_upper, eta_method, eta)

  if (!is.null(seed)) set.seed(seed)

  # 4. Calibration Phase
  if (calibrate && eta_method == "GPC") {

    # Milestone message
    message(sprintf("Calibrating eta via GPC for %s-sided %s target...", side, target))

    # Use user-supplied eta as the starting value (eta0) if provided, else default to 1.0
    eta_start <- if (!is.null(eta)) eta else 1.0

    eta_val <- .calibrate_eta_gpc(
      x = x,
      side = side,
      target = target,
      content = content,
      confidence = confidence,
      tau = tau_use,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      verbose = verbose,
      eta0 = eta_start  # <--- Pass the starting value here
    )
  } else {
    # If not calibrating, use supplied eta or default to 1.0
    if (is.null(eta)) {
      if (eta_method == "fixed") stop("Must supply 'eta' when eta_method = 'fixed'.")
      eta_val <- 1.0
    } else {
      eta_val <- eta
    }
  }
  # 5. Final Sampling Phase
  # Milestone message: Always shown
  message(sprintf("Generating final posterior samples (n = %d)...", n_mcmc))

  samples <- .run_sampler(
    x = x,
    side = side,
    target = target,
    tau = tau_use,
    tau_lower = tau_lower,
    tau_upper = tau_upper,
    eta = eta_val,
    n_mcmc = n_mcmc,
    burnin = burnin,
    verbose = verbose
  )

  # 6. Final Interval Calculation
  interval <- .compute_interval(
    samples = samples,
    eta_val = eta_val,
    side = side,
    type = type,
    content = content,
    confidence = confidence
  )

  # 7. Build S3 Object
  fit <- list(
    call = match.call(),
    data = x,
    side = side,
    type = if(side == "one") type else "centered",
    target = target,
    eta = eta_val,
    tau_used = if(side == "one") tau_use else c(tau_lower, tau_upper),
    interval = interval,
    samples = samples,
    settings = list(
      content = content,
      confidence = confidence,
      n_mcmc = n_mcmc,
      eta_method = eta_method
    )
  )
  class(fit) <- "gibbsTI"

  return(fit)
}
