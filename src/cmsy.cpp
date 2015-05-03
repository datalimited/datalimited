#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List schaefer_cmsy(NumericVector r_lim, NumericVector k_lim, double sig_r,
  NumericVector startbio, NumericVector yr, NumericVector ct, int interyr_index,
  NumericVector interbio, int reps, NumericVector finalbio) {

  int nstartbio = startbio.size();
  int nyr = yr.size();
  double ell;
  double j;
  double J;
  NumericVector ri;
  NumericVector ki;
  double r;
  double k;
  NumericVector bt(nyr+1);
  NumericVector xt(nyr+1);
  NumericMatrix out(reps, 4);
  NumericMatrix biomass_out(reps, nyr+1);

  // double prior_log_mean_minus_log2 = prior_log_mean - log(2); // for speed
  int interyr_index_minus_one = interyr_index - 1; // for speed

  ri = runif(reps, r_lim(0), r_lim(1)); // CM: changed to uniform on regular scale
  ki = exp(runif(reps, log(k_lim(0)), log(k_lim(1))));

  for (int ii=0; ii<reps; ii++) {
    r = ri(ii);
    k = ki(ii);
    ell = 0;

    for (int a=0; a<nstartbio; a++) {
      j = startbio(a);
      if (ell == 0) {
        xt = rnorm(nyr + 1, 0, sig_r);
        bt(0) = j * k * exp(xt(0)); // set biomass in first year
        for (int i=0; i<nyr; i++) {
          // double xt = 0.1;
          // calculate biomass as function of previous year's biomass plus net
          // production minus catch:
          bt(i + 1) = (bt(i) + r * bt(i) * (1 - bt(i)/k) - ct(i)) * exp(xt(i));
        }
        // Bernoulli likelihood, assign 0 or 1 to each combination of r and k
        ell = 0;
        // New posterior predictive prior on final year biomass
        // double current_bio_ratio = bt(nyr)/k;
        // double tmp = R::runif(0, 1);
        // double test = R::dlnorm(current_bio_ratio,
        //     prior_log_mean_minus_log2, prior_log_sd, 0) /
        //   R::dlnorm(exp(prior_log_mean_minus_log2),
        //       prior_log_mean_minus_log2, prior_log_sd, 0);
        // if (tmp < test &&
        //     min(bt) > 0 &&
        //     max(bt) <= k &&
        //     bt(interyr_index_minus_one)/k >= interbio(0) && // -1 because C++
        //     bt(interyr_index_minus_one)/k <= interbio(1)) { // -1 because C++
        //

        if (bt(nyr)/k >= finalbio(0) && bt(nyr)/k <= finalbio(1) && min(bt) > 0 &&
          max(bt) <= k && bt(interyr_index_minus_one)/k >=
          interbio(0) && bt[interyr_index_minus_one] / k <= interbio(1)) {

          ell = 1;
        }
        J = j;
      }
    }
    out(ii, 0) = r;
    out(ii, 1) = k;
    out(ii, 2) = ell;
    out(ii, 3) = J;

    biomass_out(ii, _) = bt;
  }
  return Rcpp::List::create(
    Rcpp::Named("theta") =
      Rcpp::DataFrame::create(
        Rcpp::Named("r")       = out(_, 0),
        Rcpp::Named("k")       = out(_, 1),
        Rcpp::Named("ell")     = out(_, 2),
        Rcpp::Named("J")       = out(_, 3)),
        Rcpp::Named("biomass") = biomass_out);
}
