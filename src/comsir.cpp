#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double comsir_effort(double a, double x, double k, double bt) {
  return x * (bt/(k*a) - 1.0);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix comsir_bounds(NumericVector x_bounds,
  NumericVector a_bounds, NumericVector k_bounds,
  NumericVector effort_bounds, int nsim, bool dampen) {

  int samples = 0;
  NumericMatrix out(nsim, 5);
  Rcpp::LogicalVector within_bounds;

  while(samples < nsim) {
    double a_ = R::runif(a_bounds(0), a_bounds(1));
    double x_ = R::runif(x_bounds(0), x_bounds(1));
    double k_ = exp(R::runif(log(k_bounds(0)), log(k_bounds(1))));
    double lower = comsir_effort(a_, x_, k_, 0.0);
    double upper = comsir_effort(a_, x_, k_, k_);

    if (dampen) {
      within_bounds =
        lower >= effort_bounds(0) & upper <= effort_bounds(1) & x_ > 0 & x_ < -a_/(a_-1);
    } else {
      within_bounds =
        lower >= effort_bounds(0) & upper <= effort_bounds(1);
    }

    if (within_bounds) {
      // if (lower >= effort_bounds(0) & upper <= effort_bounds(1)) {
      out(samples, 0) = a_;
      out(samples, 1) = x_;
      out(samples, 2) = k_;
      out(samples, 3) = lower;
      out(samples, 4) = upper;
      samples = samples + 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector check_comsir_lik(int nsim, NumericVector predbio, NumericVector
  predprop, NumericVector like) {
  for (int j=0; j<nsim; j++) {
    if ((like(j) != 0.0 && (predbio(j) <= 0.0 || predprop(j) < 0.0 )) ||
      NumericVector::is_na(like(j))) {
      like(j) = 0.0;
    }
  }
  return like;
}

// @examples
// posfun(c(0.0, 0.1, 0.2))
// [[Rcpp::export]]
NumericMatrix posfun(NumericVector x, double eps = 0.00001) {
  // returns a new value to replace the old and a penalty for the likelihood
  int x_size = x.size();
  NumericMatrix out(x_size, 2);
  double y;
  double pen;
  double newvalue;

  for (int i=0; i<x_size; i++) {
    if (x(i) >= eps) {
      out(i, 0)  = x(i);
      out(i, 1)  = 0.0;
    } else {
      pen        = 0.01 * pow(x(i) - eps, 2);
      y          = 1.0 - x(i) / eps;
      newvalue   = eps * (1.0 / (1.0 + y + y * y));
      out(i, 0)  = newvalue;
      out(i, 1)  = pen;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix effortdyn(NumericVector h, NumericVector k, NumericVector r,
  NumericVector x, NumericVector a, NumericVector yr, NumericVector ct, bool
  logistic_model) {

  int nyr = ct.size();
  int nsim = h.size();
  NumericMatrix out(nyr * nsim, 6);
  NumericVector BoverBmsy(nyr);
  double BMSY;
  int row_id = 0;

  for (int i=0; i<nsim; i++) {
    NumericVector predbio(nyr, 0.0); // first arg is length, second arg is default value
    NumericVector predprop(nyr, 0.0);
    NumericVector predcatch(nyr, 0.0);

    // initial conditions
    predbio(0) = (1.0 - h(i)) * k(i);
    predprop(0) = ct(0) / predbio(0);
    predcatch(0) = ct(0);

    for (int t=1; t<nyr; t++) { // note this starts at 1 not 0
      // Note that there is no `obs` option here: this is assuming that
      // `obs` is FALSE
      // biomass dynamics
      predbio[t] = predbio(t-1) + (r(i) * predbio(t-1) * (1 - (predbio(t-1) / k(i)))) - predcatch(t-1);
      // effort dynamics
      if (logistic_model) { // logistic model with Bt-1
        predprop(t) = predprop(t-1) * (1 + x(i) * ((predbio(t-1) / (k(i) * a(i))) - 1));
      }
      else { // linear model
        predprop(t) = predprop(t-1) + x(i) * predprop(0);
      }
      predcatch(t) = predbio(t) * predprop(t);
    }
    BMSY = k(i) / 2.0;
    BoverBmsy = predbio / BMSY;

    for (int j=0; j<nyr; j++) {
      int row_id_plus_j = row_id + j;
      out(row_id_plus_j, 0) = BoverBmsy(j);
      out(row_id_plus_j, 1) = predbio(j);
      out(row_id_plus_j, 2) = predprop(j);
      out(row_id_plus_j, 3) = predcatch(j);
      out(row_id_plus_j, 4) = ct(j);
      out(row_id_plus_j, 5) = yr(j);
    }
    row_id = row_id + nyr; // set up for next slot
  }
  return out;
}
