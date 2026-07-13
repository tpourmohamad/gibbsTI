#' Fit a Gibbs Posterior Tolerance Interval
#'
#' Constructs a one-sided or two-sided nonparametric tolerance interval from a
#' calibrated Gibbs posterior. The learning rate `eta` is chosen automatically by
#' Gibbs posterior calibration (GPC) so that the interval attains the nominal
#' frequentist coverage, or it can be fixed by the user.
#'
#' @param x Numeric vector of observations (at least 5, all finite).
#' @param side Direction of the interval: `"one"` (one-sided) or `"two"`
#'   (two-sided).
#' @param type For one-sided intervals only, `"upper"` (bound above) or `"lower"`
#'   (bound below). Ignored when `side = "two"`.
#' @param target Calibration target: `"content"` (default; the interval is
#'   defined by a population proportion) or `"quantile"` (the interval targets
#'   explicit population quantiles supplied through `tau`, `tau_lower`,
#'   `tau_upper`).
#' @param content Population proportion the interval should contain, in (0, 1)
#'   (e.g. 0.95). Used when `target = "content"`; when `target = "quantile"` it
#'   only supplies defaults for any unset `tau_*` argument.
#' @param confidence Frequentist confidence level in (0, 1) (e.g. 0.90). This is
#'   both the GPC coverage target (`1 - alpha`) and the posterior level read off
#'   as the reported limit.
#' @param tau Target quantile for a one-sided `"quantile"` interval, in (0, 1).
#'   Only used when `target = "quantile"`; ignored when `target = "content"`.
#' @param tau_lower,tau_upper Lower/upper target quantiles for a two-sided
#'   `"quantile"` interval, in (0, 1) with `tau_lower < tau_upper`. Only used
#'   when `target = "quantile"`; ignored when `target = "content"`.
#' @param eta_method How the Gibbs learning rate `eta` is chosen: `"GPC"`
#'   (default; data-driven calibration) or `"fixed"` (use the supplied `eta`).
#' @param eta Gibbs learning rate (positive). Optional starting value when
#'   `eta_method = "GPC"` (defaults to 1); required when `eta_method = "fixed"`.
#' @param max_iter Maximum number of Robbins-Monro iterations for the GPC search
#'   (default 15). Only used when `eta_method = "GPC"`.
#' @param control A list of control parameters for the MCMC samplers and GPC:
#' \itemize{
#'    \item \code{B}: Number of bootstrap replicates for calibration (default 200).
#'    \item \code{tol}: Convergence tolerance for eta (default 1e-3).
#'    \item \code{c}, \code{gamma}: Step-size parameters for Robbins-Monro (defaults 0.5, 0.75).
#'    \item \code{w}, \code{m}: Slice sampler width and max steps for one-sided intervals.
#'    \item \code{prop_sd}: Proposal standard deviation for Metropolis-Hastings in the two-sided case.
#'    \item \code{n_samps_calib}: MCMC samples used per bootstrap during calibration (default 1000).
#'    \item \code{burnin_calib}: Burn-in used per bootstrap during calibration (default 200).
#' }
#' @param n_mcmc Number of posterior draws retained for the final interval, after
#'   burn-in and thinning (default 5000).
#' @param burnin Number of initial iterations discarded as burn-in (default 1000).
#' @param thin Thinning interval: every `thin`-th post-burn-in draw is kept
#'   (default 1).
#' @param seed Optional integer random seed for reproducibility.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @return An object of class \code{"gibbsTI"}: a list with the fitted
#'   `interval` (a numeric limit for one-sided fits, or a list with `lower` and
#'   `upper` for two-sided fits), the calibrated `eta`, the retained posterior
#'   `samples`, the `data`, the quantile level(s) `tau_used`, calibration
#'   diagnostics (`calibrated_coverage`, `eta_history`), and the `settings` used.
#'   Use [print()], [summary()], and [plot()] on the result.
#' @examples
#' set.seed(1)
#' x <- rnorm(50)
#'
#' # Fast: fix the learning rate (no calibration search)
#' fit <- gibbsTI(x, side = "one", type = "upper", eta_method = "fixed",
#'                eta = 1, n_mcmc = 500, burnin = 200, verbose = FALSE)
#' print(fit)
#'
#' \donttest{
#' # Data-driven GPC calibration (slower) for a two-sided 95%/90% interval
#' fit2 <- gibbsTI(x, side = "two", content = 0.95, confidence = 0.90,
#'                 seed = 1, verbose = FALSE)
#' summary(fit2)
#' }
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
    eta_method = c("GPC", "fixed"),
    eta = NULL,
    max_iter = 15,
    control = list(),
    n_mcmc = 5000,
    burnin = 1000,
    thin = 1,
    seed = NULL,
    verbose = TRUE
) {

  # 1. Match Arguments
  side <- match.arg(side)
  type <- match.arg(type)
  target <- match.arg(target)
  eta_method <- match.arg(eta_method)

  # Default control settings
  con <- list(B = 200, tol = 1e-3, c = 0.5, gamma = 0.75, w = 0.5, m = 1000, prop_sd = 0.5,
              n_samps_calib = 1000, burnin_calib = 200)
  con[names(control)] <- control

  # 2. Logic to set tau values based on target
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
    tau_use <- if(!is.null(tau)) tau else (if(type == "upper") content else (1-content))
    if (side == "two") {
      if(is.null(tau_lower)) tau_lower <- (1 - content) / 2
      if(is.null(tau_upper)) tau_upper <- 1 - (1 - content) / 2
    }
  }

  # 2b. Validate inputs (uses the resolved tau values)
  .validate_inputs(x, content, confidence, target, tau_use, tau_lower, tau_upper,
                   eta_method, eta, max_iter, n_mcmc, burnin, thin)

  if (!is.null(seed)) set.seed(seed)

  # Variables to store calibration diagnostics
  final_coverage <- NULL
  eta_history <- NULL

  # 3. Calibration Phase
  if (eta_method == "GPC") {
    if(verbose) message(sprintf("Calibrating eta via GPC for %s-sided %s target...", side, target))

    eta_start <- if (!is.null(eta)) eta else 1.0

    # Capture the full list from C++
    calib_res <- .calibrate_eta_gpc(
      x = x,
      side = side,
      target = target,
      content = content,
      confidence = confidence,
      tau = tau_use,
      tau_lower = tau_lower,
      tau_upper = tau_upper,
      verbose = verbose,
      eta0 = eta_start,
      max_iter = max_iter,
      B = con$B,
      tol = con$tol,
      c = con$c,
      gamma = con$gamma,
      prop_sd = con$prop_sd,
      w = con$w,
      m = con$m,
      n_samps_boot = con$n_samps_calib,
      burn_in_boot = con$burnin_calib
    )

    # Extract values from the C++ returned list
    eta_val        <- calib_res$final_eta
    final_coverage <- calib_res$final_coverage
    eta_history    <- calib_res$eta_history

  } else {
    # eta_method == "fixed": use the user-supplied learning rate
    # (.validate_inputs() has already checked that a positive eta was provided).
    eta_val <- eta
  }

  # 4. Final Sampling Phase
  if(verbose) message(sprintf("Generating final posterior samples (n = %d)...", n_mcmc))

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
    thin = thin,
    w = con$w,
    m = con$m,
    prop_sd = con$prop_sd,
    verbose = verbose
  )

  # 5. Final Interval Calculation
  interval <- .compute_interval(
    samples = samples,
    eta_val = eta_val,
    side = side,
    type = type,
    content = content,
    confidence = confidence
  )

  # 6. Build S3 Object
  fit <- list(
    call = match.call(),
    data = x,
    side = side,
    type = if(side == "one") type else "centered",
    target = target,
    eta = eta_val,
    calibrated_coverage = final_coverage,
    eta_history = eta_history,
    tau_used = if(side == "one") tau_use else c(tau_lower, tau_upper),
    interval = interval,
    samples = samples,
    settings = list(
      content = content,
      confidence = confidence,
      n_mcmc = n_mcmc,
      eta_method = eta_method,
      control = con
    )
  )
  class(fit) <- "gibbsTI"

  return(fit)
}
