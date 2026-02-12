#include <Rcpp.h>
#include <iomanip>

using namespace Rcpp;

// Utility: Check loss
double rho_check_cpp(NumericVector r, double tau) {
  int n = r.size();
  double loss = 0.0;
  for(int i = 0; i < n; i++) {
    loss += r[i] * (tau - (r[i] < 0 ? 1.0 : 0.0));
  }
  return loss;
}

double logpost_q_cpp(double q, NumericVector y, double tau, double eta) {
  NumericVector r = y - q;
  return -eta * rho_check_cpp(r, tau);
}

// [[Rcpp::export]]
NumericVector slice_sample_1d_cpp(double x0, NumericVector y, double tau, double eta,
                                  double w = 0.1, int m = 100, int n_samps = 4000,
                                  double lower = R_NegInf, double upper = R_PosInf) {
  NumericVector x(n_samps);
  double x_cur = x0;
  for(int i = 0; i < n_samps; i++) {
    double logy = logpost_q_cpp(x_cur, y, tau, eta) - R::rexp(1.0);
    double u = R::runif(0.0, w);
    double L = x_cur - u;
    double R = x_cur + (w - u);
    int j = floor(m * R::runif(0.0, 1.0));
    int k = (m - 1) - j;
    while(j > 0 && L > lower && logpost_q_cpp(L, y, tau, eta) > logy) { L -= w; j--; }
    while(k > 0 && R < upper && logpost_q_cpp(R, y, tau, eta) > logy) { R += w; k--; }
    L = std::max(L, lower);
    R = std::min(R, upper);
    while(true) {
      double x_prop = R::runif(L, R);
      double logf_prop = logpost_q_cpp(x_prop, y, tau, eta);
      if(std::isfinite(logf_prop) && logf_prop >= logy) { x_cur = x_prop; break; }
      if(x_prop < x_cur) L = x_prop; else R = x_prop;
    }
    x[i] = x_cur;
  }
  return x;
}

// [[Rcpp::export]]
double bootstrap_coverage_cpp(double tau, double alpha, int B, double eta, double theta_hat,
                              NumericMatrix boots, int n_samps_boot, int burn_in_boot,
                              double w, int m) {
  Function quantile("quantile");
  int n_total_mcmc = n_samps_boot + burn_in_boot;
  NumericVector coverage_indicators(B);
  double lower_prob = alpha / 2.0;
  double upper_prob = 1.0 - alpha / 2.0;

  for(int b = 0; b < B; b++) {
    NumericVector boot_data = boots(_, b);

    // Sample with calibration-specific settings
    NumericVector boot_draws = slice_sample_1d_cpp(Rcpp::median(boot_data), boot_data, tau, eta, w, m, n_total_mcmc);

    // Manual Burn-in removal (as Rcpp doesn't have a slice operator like R)
    NumericVector boot_posterior(n_samps_boot);
    for(int i = 0; i < n_samps_boot; i++) {
      boot_posterior[i] = boot_draws[i + burn_in_boot];
    }

    NumericVector probs = NumericVector::create(lower_prob, upper_prob);
    NumericVector cred_interval = quantile(boot_posterior, Named("probs") = probs);
    coverage_indicators[b] = (theta_hat >= cred_interval[0] && theta_hat <= cred_interval[1]) ? 1.0 : 0.0;
  }
  return mean(coverage_indicators);
}

// [[Rcpp::export]]
List calibrate_eta_cpp(NumericVector y, double tau, double alpha = 0.1, double eta0 = 1.0,
                       int B = 200, int max_iter = 15, double tol = 1e-3,
                       double c = 0.5, double gamma = 0.75,
                       int n_samps_boot = 1000, int burn_in_boot = 200,
                       double w = 0.1, int m = 100, bool verbose = true) {
  int n = y.size();

  // Initial estimate (using n=4000 for a stable "truth")
  NumericVector initial_posterior = slice_sample_1d_cpp(Rcpp::median(y), y, tau, 1.0, w, m, 4000);
  double theta_hat = Rcpp::median(initial_posterior);

  NumericMatrix boots(n, B);
  for(int b = 0; b < B; b++) boots(_, b) = sample(y, n, true);

  double eta_current = eta0;
  NumericVector eta_history(max_iter + 1);
  eta_history[0] = eta_current;

  int s;
  for(s = 1; s <= max_iter; s++) {
    double kappa_s = c / std::pow(1.0 + s, gamma);
    double est_coverage = bootstrap_coverage_cpp(tau, alpha, B, eta_current, theta_hat,
                                                 boots, n_samps_boot, burn_in_boot, w, m);

    if(verbose) {
      Rcpp::Rcout << "Iteration " << s << ": eta = " << std::fixed << std::setprecision(4)
                  << eta_current << " | coverage = " << est_coverage << std::endl;
    }

    double eta_prev = eta_current;
    eta_current = eta_prev + kappa_s * (est_coverage - (1.0 - alpha));

    if(eta_current <= 0) eta_current = 1e-4;
    eta_history[s] = eta_current;

    if(std::abs(eta_current - eta_prev) < tol) {
      if(verbose) Rcpp::Rcout << "Converged early at iteration " << s << std::endl;
      break;
    }
  }

  return List::create(_["final_eta"] = eta_current,
                      _["eta_history"] = eta_history[Range(0, (s < max_iter ? s : (max_iter)))]);
}
