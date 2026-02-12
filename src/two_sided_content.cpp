#include <Rcpp.h>
#include <iomanip>

using namespace Rcpp;

// ----------------------------
// Utility: Check loss
// ----------------------------
static double rho_check_internal(NumericVector r, double tau) {
  int n = r.size();
  double loss = 0.0;

  for(int i = 0; i < n; i++) {
    loss += r[i] * (tau - (r[i] < 0 ? 1.0 : 0.0));
  }

  return loss;
}

// ----------------------------
// Log posterior
// ----------------------------
double logpost_joint_content(double t1,
                             double t2,
                             NumericVector y,
                             double tl,
                             double tu,
                             double eta) {

  double ql = t1;
  double qu = t1 + exp(t2);

  return -eta *
    (rho_check_internal(y - ql, tl) +
    rho_check_internal(y - qu, tu)) +
    t2;
}

// ----------------------------
// MCMC sampler (Metropolis-Hastings)
// ----------------------------

// [[Rcpp::export]]
NumericMatrix sample_joint_cpp(NumericVector y,
                               double tau_lower,
                               double tau_upper,
                               double eta,
                               int n_samps = 4000,
                               int burn_in = 1000,
                               double prop_sd = 0.5) {

  int total = n_samps + burn_in;
  NumericMatrix out(n_samps, 2);

  double t1 = Rcpp::median(y);
  double t2 = 0.1;
  double logp = logpost_joint_content(t1, t2, y, tau_lower, tau_upper, eta);

  for(int i = 0; i < total; i++) {
    // Random walk proposal using prop_sd
    double t1p = R::rnorm(t1, prop_sd);
    double t2p = R::rnorm(t2, prop_sd);

    double logpp = logpost_joint_content(t1p, t2p, y, tau_lower, tau_upper, eta);

    if(log(R::runif(0,1)) < logpp - logp) {
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

// ----------------------------
// Bootstrap coverage
// ----------------------------

// [[Rcpp::export]]
double bootstrap_coverage_content_cpp(double tl,
                                      double tu,
                                      double alpha,
                                      int B,
                                      double eta,
                                      NumericMatrix boots,
                                      NumericVector y_orig,
                                      int n_samps_boot,  // Added to control speed
                                      int burn_in_boot,  // Added to control speed
                                      double prop_sd) {   // Added to match sampler

  // Access the package namespace for the interval calculation
  Environment pkg = Environment::namespace_env("gibbsTI");
  Function compute_interval = pkg["compute_two_sided_interval"];

  double P = tu - tl;
  double successes = 0.0;

  for(int b = 0; b < B; b++) {
    // Run the sampler on the bootstrap data
    NumericMatrix post = sample_joint_cpp(boots(_, b), tl, tu, eta, n_samps_boot, burn_in_boot, prop_sd);

    // Calculate the interval from posterior samples
    List interval = compute_interval(post(_,0), post(_,1), 1.0 - alpha);

    double L = interval["lower"];
    double U = interval["upper"];

    // Check how many original data points fall within this bootstrap interval
    int cnt = 0;
    for(int i = 0; i < y_orig.size(); i++) {
      if(y_orig[i] >= L && y_orig[i] <= U)
        cnt++;
    }

    if(((double)cnt / y_orig.size()) >= P)
      successes += 1.0;
  }

  return successes / B;
}

// ----------------------------
// Robbins–Monro calibration
// ----------------------------

// [[Rcpp::export]]
List calibrate_eta_joint_content_cpp(NumericVector y,
                                     double tau_lower,
                                     double tau_upper,
                                     double alpha = 0.1,
                                     double eta0 = 1.0,
                                     int B = 200,
                                     int max_iter = 15,
                                     double tol = 1e-3, // Added tol to match R wrapper
                                     double c = 0.5,
                                     double gamma = 0.75,
                                     double prop_sd = 0.5,
                                     int n_samps_boot = 1000, // <--- ADD THIS
                                     int burn_in_boot = 200,  // <--- ADD THIS
                                     bool verbose = true) {

  int n = y.size();
  NumericMatrix boots(n, B);

  for(int b = 0; b < B; b++) {
    boots(_, b) = sample(y, n, true);
  }

  double eta = eta0;
  double eta_prev = eta;

  for(int s = 1; s <= max_iter; s++) {
    // Pass the actual n_samps_boot and burn_in_boot variables here!
    double cover = bootstrap_coverage_content_cpp(tau_lower, tau_upper, alpha, B, eta,
                                                  boots, y, n_samps_boot, burn_in_boot, prop_sd);

    if(verbose) {
      Rcpp::Rcout << "Iteration " << s << ": eta = " << std::fixed << std::setprecision(4)
                  << eta << " | coverage = " << cover << std::endl;
    }

    eta_prev = eta;
    double step = c / pow(1.0 + s, gamma);
    eta += step * (cover - (1.0 - alpha));

    if(eta <= 0) eta = 1e-4;

    // Optional: Add early convergence check to match your one-sided logic
    if(std::abs(eta - eta_prev) < tol) {
      if(verbose) Rcpp::Rcout << "Converged early at iteration " << s << std::endl;
      break;
    }
  }

  return List::create(_["final_eta"] = eta);
}
