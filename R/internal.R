.validate_inputs <- function(x, content, confidence, target, tau_use, tau_lower, tau_upper, eta_method, eta) {

  # 1. Data Checks (x)
  if (!is.numeric(x))
    stop("x must be a numeric vector.")

  if (length(x) < 5)
    stop("x must contain at least 5 observations for meaningful intervals.")

  if (any(!is.finite(x)))
    stop("x contains non-finite values (NA, NaN, or Inf).")

  # 2. Probability Checks (content & confidence)
  if (!is.numeric(content) || length(content) != 1 || content <= 0 || content >= 1)
    stop("content must be a single number in (0, 1).")

  if (!is.numeric(confidence) || length(confidence) != 1 || confidence <= 0 || confidence >= 1)
    stop("confidence must be a single number in (0, 1).")

  # 3. Target / Tau Logic
  if (target == "quantile") {
    # Check for one-sided tau_use or two-sided tau pair
    if (is.null(tau_use) && (is.null(tau_lower) || is.null(tau_upper))) {
      stop("Either tau (one-sided) or tau_lower/upper (two-sided) must be supplied when target = 'quantile'.")
    }

    # Validate range of supplied taus
    all_taus <- c(tau_use, tau_lower, tau_upper)
    if (any(all_taus <= 0 | all_taus >= 1))
      stop("All tau values must lie in the interval (0, 1).")

    # Logic check: Lower must be less than Upper
    if (!is.null(tau_lower) && !is.null(tau_upper)) {
      if (tau_lower >= tau_upper) {
        stop("tau_lower must be strictly less than tau_upper.")
      }
    }
  }

  # 4. Learning Rate Logic (eta)
  if (eta_method == "fixed" && is.null(eta)) {
    stop("When eta_method is 'fixed', a numeric 'eta' value must be provided.")
  }

  if (!is.null(eta)) {
    if (!is.numeric(eta) || length(eta) != 1 || eta <= 0) {
      stop("eta must be a single positive numeric value.")
    }
  }

  invisible(TRUE)
}
