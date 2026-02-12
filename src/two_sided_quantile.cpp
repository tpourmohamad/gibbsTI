#include <Rcpp.h>
#include <iomanip>

using namespace Rcpp;

// =============================================================================
// UTILITIES
// =============================================================================

static double rho_check_internal(NumericVector r, double tau) {
  int n = r.size();
  double loss = 0.0;

  for(int i = 0; i < n; i++) {
    loss += r[i] * (tau - (r[i] < 0 ? 1.0 : 0.0));
  }

  return loss;
}

double logpost_joint_quantile(double theta1,
                              double theta2,
                              NumericVector y,
                              double tau_lower,
                              double tau_upper,
                              double eta) {

  double q_lower = theta1;
  double q_upper = theta1 + exp(theta2);

  double joint_loss =
    rho_check_internal(y - q_lower, tau_lower) +
    rho_check_internal(y - q_upper, tau_upper);

  return -eta * joint_loss + theta2;
}

// =============================================================================
// 1. 2D SAMPLER (Metropolis-Hastings)
// =============================================================================

// [[Rcpp::export]]
NumericMatrix sample_joint_quantile_cpp(NumericVector y,
                                        double tau_lower,
                                        double tau_upper,
                                        double eta,
                                        int n_samps = 4000,
                                        int burn_in = 1000,
                                        double prop_sd = 0.5) {

  int total_samps = n_samps + burn_in;
  NumericMatrix out(n_samps, 2);

  double t1 = Rcpp::median(y);
  double t2 = 0.1;

  double logp = logpost_joint_quantile(t1, t2, y, tau_lower, tau_upper, eta);

  for(int i = 0; i < total_samps; i++) {
    double t1p = R::rnorm(t1, prop_sd);
    double t2p = R::rnorm(t2, prop_sd);

    double logpp = logpost_joint_quantile(t1p, t2p, y, tau_lower, tau_upper, eta);

    if(std::isfinite(logpp - logp) && log(R::runif(0,1)) < logpp - logp) {
      t1 = t1p;
      t2 = t2p;
      logp = logpp;
    }

    if(i >= burn_in) {
      out(i - burn_in, 0) = t1;
      out(i - burn_in, 1) = t1 + exp(t2);
    }
  }

  return out;
}

// =============================================================================
// 2. QUANTILE-BASED BOOTSTRAP COVERAGE
// =============================================================================

// [[Rcpp::export]]
double bootstrap_coverage_quantile_cpp(double tl,
                                       double tu,
                                       double alpha,
                                       int B,
                                       double eta,
                                       NumericVector theta_hat_vec,
                                       NumericMatrix boots,
                                       int n_samps_boot,
                                       int burn_in_boot,
                                       double prop_sd) {

  Environment pkg = Environment::namespace_env("gibbsTI");
  Function compute_interval = pkg["compute_two_sided_interval"];

  NumericVector coverage_indicators(B);

  for(int b = 0; b < B; b++) {
    NumericMatrix post = sample_joint_quantile_cpp(boots(_, b), tl, tu, eta, n_samps_boot, burn_in_boot, prop_sd);

    List interval = compute_interval(
      _["q_lower"]    = post(_, 0),
      _["q_upper"]    = post(_, 1),
      _["conf_level"] = 1.0 - alpha
    );

    double L = interval["lower"];
    double U = interval["upper"];

    bool covered = (theta_hat_vec[0] >= L && theta_hat_vec[1] <= U);
    coverage_indicators[b] = covered ? 1.0 : 0.0;
  }

  return mean(coverage_indicators);
}

// =============================================================================
// 3. CALIBRATION FOR QUANTILE MODEL
// =============================================================================

// [[Rcpp::export]]
List calibrate_eta_joint_quantile_cpp(NumericVector y,
                                      double tau_lower,
                                      double tau_upper,
                                      double alpha = 0.1,
                                      double eta0 = 1.0,
                                      int B = 200,
                                      int max_iter = 15,
                                      double tol = 1e-3,
                                      double c = 0.5,
                                      double gamma = 0.75,
                                      double prop_sd = 0.5,
                                      int n_samps_boot = 1000,
                                      int burn_in_boot = 200,
                                      bool verbose = true) {

  int n = y.size();

  // Initial point estimates
  NumericMatrix init_post = sample_joint_quantile_cpp(y, tau_lower, tau_upper, 1.0, 2000, 500, prop_sd);

  NumericVector col0 = init_post(_, 0);
  NumericVector col1 = init_post(_, 1);

  NumericVector theta_hat_vec = NumericVector::create(
    Rcpp::median(col0),
    Rcpp::median(col1)
  );

  NumericMatrix boots(n, B);
  for(int b = 0; b < B; b++) {
    boots(_, b) = sample(y, n, true);
  }

  double eta = eta0;
  double eta_prev = eta;
  NumericVector hist(max_iter + 1);
  hist[0] = eta;

  int s;
  for(s = 1; s <= max_iter; s++) {
    double kappa = c / std::pow(1.0 + s, gamma);

    double cover = bootstrap_coverage_quantile_cpp(tau_lower, tau_upper, alpha, B, eta,
                                                   theta_hat_vec, boots, n_samps_boot, burn_in_boot, prop_sd);

    if(verbose) {
      Rcpp::Rcout << "Iteration " << s << ": eta = " << std::fixed << std::setprecision(4)
                  << eta << " | coverage = " << cover << std::endl;
    }

    eta_prev = eta;
    eta += kappa * (cover - (1.0 - alpha));

    if(eta <= 0) eta = 1e-4;
    hist[s] = eta;

    if(std::abs(eta - eta_prev) < tol) {
      if(verbose) Rcpp::Rcout << "Converged early at iteration " << s << std::endl;
      break;
    }
  }

  return List::create(
    _["final_eta"]   = eta,
    _["eta_history"] = hist[Range(0, (s < max_iter ? s : max_iter))],
                               _["theta_hat_vec"] = theta_hat_vec
  );
}
