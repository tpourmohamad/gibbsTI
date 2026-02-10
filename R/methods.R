#' @export
print.gibbsTI <- function(x, ...) {
  int_val <- unlist(x$interval)

  if (x$side == "one") {
    label <- if (!is.null(x$type)) tools::toTitleCase(x$type) else "One-sided"
    cat(label, "Gibbs TI:", round(int_val, 4), "\n")
  } else {
    cat("Two-sided Gibbs TI: [", round(int_val[1], 4), ", ",
        round(int_val[2], 4), "]\n", sep="")
  }
  cat("Call summary() for full details.\n")
  invisible(x)
}

#' @export
summary.gibbsTI <- function(object, ...) {
  cat("\n==========================================")
  cat("\n  Gibbs Posterior Tolerance Interval")
  cat("\n==========================================\n")

  # Calculate effective content
  if (object$side == "one") {
    # For one-sided, it's either P(X < U) or P(X > L)
    eff_content <- if (object$type == "upper") object$tau_used else (1 - object$tau_used)
  } else {
    # For two-sided, it's the mass between the two taus
    eff_content <- object$tau_used[2] - object$tau_used[1]
  }

  cat(sprintf("%-12s: %s\n", "Side", toupper(object$side)))
  cat(sprintf("%-12s: %s\n", "Target", tools::toTitleCase(object$target)))
  cat(sprintf("%-12s: %g%%\n", "Content", eff_content * 100))
  cat(sprintf("%-12s: %g%%\n", "Confidence", object$settings$confidence * 100))

  cat("------------------------------------------\n")
  if (object$side == "one") {
    label <- if (object$type == "upper") "Upper TI" else "Lower TI"
    cat(sprintf("%-12s: %f\n", label, unlist(object$interval)))
  } else {
    int_val <- unlist(object$interval)
    cat(sprintf("%-12s: [%.4f, %.4f]\n", "Interval", int_val[1], int_val[2]))
  }

  cat(sprintf("%-12s: %.6f\n", "Final Eta", object$eta))
  cat("==========================================\n")
  invisible(object)
}

#' @export
plot.gibbsTI <- function(x, ...) {
  dat <- x$data
  int_val <- unlist(x$interval)

  # Calculate effective content for the title
  if (x$side == "one") {
    eff_content <- if (x$type == "upper") x$tau_used else (1 - x$tau_used)
  } else {
    eff_content <- x$tau_used[2] - x$tau_used[1]
  }

  hist(dat,
       main = paste0(round(eff_content*100), "% Content / ",
                     round(x$settings$confidence*100), "% Confidence TI"),
       xlab = "Value", col = "gray90", border = "white")

  abline(v = int_val, col = "firebrick", lwd = 2, lty = 2)

  leg_text <- if(x$side == "one") paste(tools::toTitleCase(x$type), "Limit") else "Tolerance Interval"
  legend("topright", legend = leg_text, col = "firebrick", lwd = 2, lty = 2, bty = "n")
}
