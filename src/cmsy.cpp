#include <Rcpp.h>
using namespace Rcpp;

//' @param r_lim Numeric vector of lower and upper intrinsic population growth
//'   rates
//' @param k_lim Numeric vector of lower and upper carrying capacities
//' @param startbio A vector of possible starting biomasses to loop over
//'
//' @useDynLib datalimited
//' @importFrom Rcpp sourceCpp
//' @export
//' @rdname cmsy
// [[Rcpp::export]]
DataFrame schaefer_cmsy(NumericVector r_lim, NumericVector k_lim, double sig_r,
  NumericVector startbio, NumericVector yr, NumericVector ct, int interyr_index,
  double prior_log_mean, double prior_log_sd,
  NumericVector interbio, int reps) {

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
  NumericMatrix out(reps, 4);

  double prior_log_mean_minus_log2 = prior_log_mean - log(2); // for speed
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
        bt(0) = j * k * exp(R::rnorm(0, sig_r)); // set biomass in first year
        for (int i=0; i<nyr; i++) {
          double xt = R::rnorm(0, sig_r);
          // double xt = 0.1;
          // calculate biomass as function of previous year's biomass plus net
          // production minus catch:
          bt(i + 1) = (bt(i) + r * bt(i) * (1 - bt(i)/k) - ct(i)) * exp(xt);
        }
        // Bernoulli likelihood, assign 0 or 1 to each combination of r and k
        ell = 0;
        // New posterior predictive prior on final year biomass
        double current_bio_ratio = bt(nyr)/k;
        double tmp = R::runif(0, 1);
        double test = R::dlnorm(current_bio_ratio,
          prior_log_mean_minus_log2, prior_log_sd, 0) /
        R::dlnorm(exp(prior_log_mean_minus_log2),
          prior_log_mean_minus_log2, prior_log_sd, 0);
        if (tmp < test &&
            min(bt) > 0 &&
            max(bt) <= k &&
            bt(interyr_index_minus_one)/k >= interbio(0) && // -1 because C++
            bt(interyr_index_minus_one)/k <= interbio(1)) { // -1 because C++
              ell = 1;
        }
        J = j;
      }
    }
    out(ii, 0) = r;
    out(ii, 1) = k;
    out(ii, 2) = ell;
    out(ii, 3) = J;
  }
  return DataFrame::create(
    Named("r")       = out(_, 0),
    Named("k")       = out(_, 1),
    Named("ell")     = out(_, 2),
    Named("J")       = out(_, 3));
}

// @return A matrix: each column is an iteration of the algorithm and each row
//   is a year of biomass
// [[Rcpp::export]]
NumericMatrix get_cmsy_biomass(NumericVector r, NumericVector k, NumericVector j,
  double sigR, int nyr, NumericVector ct) {
  NumericMatrix BT(nyr + 1, r.size());
  NumericVector bt(nyr + 1);
  NumericVector xt(nyr);

  for (int v=0; v<r.size(); v++) {
    xt = rnorm(nyr, 0, sigR);
    bt(0) = j(v) * k(v) * exp(R::rnorm(0, sigR)); // set biomass in first year

  for (int i=0; i<nyr; i++) { // for all years in the time series
      // calculate biomass as function of previous year's biomass plus
      // net production minus catch:
      bt(i+1) = (bt(i) + r(v) * bt(i) * (1 - bt(i)/k(v)) - ct(i)) * exp(xt(i));
    }
    BT(_, v) = bt; // exclude the initial year
  }
  return BT;
}

// @examples
// schaefer_cmsy(r_lim = c(0.1, 0.3), k_lim = c(80, 100), sig_r = 0.1,
//   startbio = seq(0.2, 0.6, by = 0.05), yr = 1:10, ct = rnorm(10),
//   prior_log_mean = 0.3, prior_log_sd = 0.1, interyr_index = 2,
//   interbio = c(0, 1), reps = 10L)
