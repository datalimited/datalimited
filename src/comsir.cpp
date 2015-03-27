#include <Rcpp.h>
using namespace Rcpp;

// @examples
//  ct <- c(512, 865, 511, 829, 935, 1390, 1260, 2240, 3370, 2770, 3050,
//    3290, 4540, 3300, 3500, 3190, 2880, 3490, 5670, 6310, 9550, 8700,
//    9130, 9160, 8490, 6400, 4420, 3680, 3190, 3960, 3290, 4220, 4220)
//
// out <- comsir_priors(ct = ct,
//    k = 800, r = 0.6, x = 0.5, a = 0.8, start_r = c(0.2, 1),
//    mink = max(ct), maxk = max(ct) * 100, nsim = 1e5))
// [[Rcpp::export]]
DataFrame comsir_priors(NumericVector ct, double k, double r, double x,
  double a, NumericVector start_r,
  double mink, double maxk,
  bool logk = true,
  double cv = 0.4, bool norm_k = false, bool norm_r = false,
  bool norm_a = false, bool norm_x = false,
  bool logistic_model = true, bool obs = false, int nsim = 2000L) {

  NumericVector k_vec(nsim);
  NumericVector r_vec(nsim);
  NumericVector a_vec(nsim);
  NumericVector x_vec(nsim);
  NumericVector h(nsim);
  NumericVector z(nsim);
  NumericVector n1(nsim);
  NumericVector like(nsim);
  NumericVector predbio(nsim);
  NumericVector predcatch(nsim);
  NumericVector predprop(nsim);
  NumericVector inipredprop(nsim);

  if (logk) {
    k_vec = exp(runif(nsim, log(mink), log(maxk)));
  } else {
    k_vec = runif(nsim, mink, maxk); // prior on arithmetic scale for K
  }

  if (!logk && norm_k) { // normal prior on the arithmetic scale
    k_vec = rnorm(nsim, k, cv * k);
  }

  if (norm_r) {
    r_vec = rnorm(nsim, r, cv * r);
  } else {
    r_vec = runif(nsim, start_r(0), start_r(1)); // prior for r
  }

  if (norm_x) {
    x_vec = rnorm(nsim, x, cv * x);
  } else {
    x_vec = runif(nsim, 0.000001, 1.0); // prior for x
  }

  if (norm_a) {
    a_vec = rnorm(nsim, a, cv * a);
  } else {
    a_vec = runif(nsim, 0.0, 1.0); // prior for a
  }

  h = rep(0.0, nsim); // h=0, start from unexploited state
  z = rep(1.0, nsim); // z=1 is the Schaeffer model
  n1 = (1.0 - h) * k_vec; // initial population size

  // like=1 if the model does not crashes, like=0 if the model crashes
  like = rep(1.0, nsim);

  // Set up the numbers and project to the start of the first "real" year
  // (i.e. n1 to Nm), all simulations at once
  // initial conditions
  predbio = n1;
  predcatch = rep(ct(0), nsim); // the first year of catches in assumed known
  predprop = predcatch / predbio;
  inipredprop = predprop;

  for (int i=0; i<(ct.size()-1); i++) { // TODO check this for errors

    // effort dynamics
    if (logistic_model) { // logistic model with Bt-1
      predprop = predprop * (1.0 + x_vec * ((predbio / (a_vec * k_vec) - 1.0)));
    } else { //linear model
      predprop = predprop + x_vec * inipredprop;
    }

    // biomass dynamics
    if (obs) {
      predbio = predbio + (r_vec * predbio * ( 1.0 - (predbio / k_vec ))) - ct(i);
    } else {
      predbio = predbio + (r_vec * predbio * ( 1.0 - (predbio / k_vec ))) - predcatch;
    }

    predcatch = predbio * predprop;

    for (int j=0; j<nsim; j++) {
      if ((like(j) != 0.0 && (predbio(j) <= 0.0 || predprop(j) < 0.0 )) ||
          NumericVector::is_na(like(j))) {
        like(j) = 0.0;
      }
    }
  }

  // in some of the years the population crashes, we know this is impossible,
  // because the population is still there (because of the catches)

  return DataFrame::create(
      Named("n1")      = n1,
      Named("k")       = k_vec,
      Named("r")       = r_vec,
      Named("z")       = z,
      Named("a")       = a_vec,
      Named("x")       = x_vec,
      Named("h")       = h,
      Named("biomass") = predbio,
      Named("prop")    = predprop,
      Named("like")    = like);
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

// @examples
//  ct <- c(512, 865, 511, 829, 935, 1390, 1260, 2240, 3370, 2770, 3050,
//    3290, 4540, 3300, 3500, 3190, 2880, 3490, 5670, 6310, 9550, 8700,
//    9130, 9160, 8490, 6400, 4420, 3680, 3190, 3960, 3290, 4220, 4220)
//
// o <- comsir_priors(ct = ct,
//    k = 800, r = 0.6, x = 0.5, a = 0.8, start_r = c(0.2, 1),
//    mink = max(ct), maxk = max(ct) * 100, nsim = 1e3)
// o <- o[o$like != 0, ]
// out <- with(o, comsir_est(n1 = n1, k = k, r = r, a = a, x = x, h = h, z = z,
//    like = like, ct))
// [[Rcpp::export]]
DataFrame comsir_est(NumericVector n1,
    NumericVector k,
    NumericVector r,
    NumericVector a,
    NumericVector x,
    NumericVector h,
    NumericVector z,
    NumericVector like,
    NumericVector ct,
    double cv = 0.4,
    bool logistic_model = false,
    bool normal_like = true) {

  // normal_like=true for normal likelihood
  // normal_like=false for lognormal likelihood
  //
  // The log normal distribution has density
  //  f(x) = 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))
  //     where mu and sigma are the mean and standard deviation of the
  //     logarithm. The mean is E(X) = exp(mu + 1/2 sigma^2), and the
  //     variance Var(X) = exp(2*mu + sigma^2)*(exp(sigma^2) - 1) and hence
  //     the coefficient of variation is sqrt(exp(sigma^2) - 1) which is
  //     approximately sigma when that is small (e.g., sigma < 1/2).
  //    x<-seq(0.001,100,0.001)
  //    y<-dlnorm(log(x),meanlog=log(1),sdlog=1,log=FALSE)
  //     dlnorm(1) == dnorm(0)
  //    z<-dnorm(log(x),mean=0,sdlog=1,log=FALSE)
  //    plot(log(x),y,type="l")
  //    plot(log(x),z,type="l")

  int nsim = k.size();
  NumericVector predbio(nsim);
  NumericVector predcatch(nsim);
  NumericVector predprop(nsim);
  NumericVector inipredprop(nsim);
  NumericVector loglike(nsim);
  NumericVector pen1(nsim);
  NumericVector pen2(nsim);
  NumericMatrix temp1(nsim, 2);
  NumericMatrix temp2(nsim, 2);
  NumericVector cum_loglike(nsim);
  NumericVector my_sd(nsim);
  double Totlike;

  // the parameters come from the resampling part
  // sigma = rnorm(nsim,0.001,0.0001)//  normal informative prior for sigma
  // Set up the numbers and project to the start of the first "real" year
  // (i.e. n1 to Nm), all simulations at once
  // initial conditions
  loglike    = rep(0.0, nsim);
  pen1       = rep(0.0, nsim);
  pen2       = rep(0.0, nsim);
  cum_loglike = rep(0.0, nsim);

  for (int i=0; i<(ct.size()-i); i++) { // TODO check this for errors
    if (i==0) {
      predbio = n1;
      predcatch = rep(ct(0), nsim); // first catch is assumed known -- BIG ASSUMPTION ---
      predprop = predcatch / predbio;
      inipredprop = predprop;
    } else {
      // effort dynamics
      if (logistic_model) { // logistic model with Bt-1
        temp1 = posfun(predprop * (1 + x * ((predbio / (a * k) - 1))));
      } else { // linear model
        temp1 = posfun(predprop + x * inipredprop);
      }
      for (int j=0; j<nsim; j++) {
        if (temp1(j, 0) > 0.99) {
          predprop(j) = 0.99;
        } else {
          predprop(j) = temp1(j, 0);
        }
        pen1(j) = pen1(j) + temp1(j, 1);
      }

      // biomass dynamics
      temp2 = posfun(predbio + (r * predbio * (1 - (predbio / k ))) - predcatch);
      for (int j=0; j<nsim; j++) {
        predbio(j) = temp2(j, 0);
        pen2(j) = pen2(j) + temp2(j, 1);
        predcatch(j) = predbio(j) * predprop(j);
      }
    }
    //  assumption of cv=0.4 in Vasconcellos and Cochrane
    if (normal_like) {
      my_sd = cv * predcatch;
      for (int j=0; j<nsim; j++) {
        loglike(j) = R::dnorm(ct(i), predcatch(j), my_sd(j), 1); // normal Log likelihood
        cum_loglike(j) = cum_loglike(j) + loglike(j); // cummulative normal loglikelihood
      }
    } else {
      for (int j=0; j<nsim; j++) {
        loglike(j) = R::dlnorm(ct(i), log(predcatch(j)), cv, 1);
        cum_loglike(j) = cum_loglike(j) + loglike(j);
      }
    }
  }
  loglike = rep(0.0, nsim); //  set vector to 0
  loglike = cum_loglike - pen1 - pen2; // the penalties turn out too big
  Totlike = sum(loglike);
  like    = exp(loglike) / sum(exp(loglike)); // re-normalize so it will add up to 1

  // in some of the years the population crashes, we know this is impossible,
  // because the population is still there (because of the catches

  return DataFrame::create(
      Named("n1")       = n1,
      Named("k")        = k,
      Named("r")        = r,
      Named("z")        = z,
      Named("a")        = a,
      Named("x")        = x,
      Named("h")        = h,
      Named("B")        = predbio,    // biomass in the latest year
      Named("prop")     = predprop,   // harvest rate in the latest year
      Named("loglike")  = cum_loglike, // log likelihood used for DIC
      Named("like")     = like);      // normalized likelihood for SIR
}

// [[Rcpp::export]]
NumericMatrix effortdyn(NumericVector h, NumericVector k, NumericVector r,
    NumericVector x, NumericVector a, NumericVector yrs, NumericVector ct, bool
    logistic_model) {

  int nyrs = ct.size();
  int nsim = h.size();
  NumericMatrix out(nyrs * nsim, 6);
  NumericVector BoverBmsy(nyrs);
  double BMSY;
  int row_id = 0;

  for (int i=0; i<nsim; i++) {
    NumericVector predbio(nyrs, 0.0); // first arg is length, second arg is default value
    NumericVector predprop(nyrs, 0.0);
    NumericVector predcatch(nyrs, 0.0);

    // initial conditions
    predbio(0) = (1.0 - h(i)) * k(i);
    predprop(0) = ct(0) / predbio(0);
    predcatch(0) = ct(0);

    for (int t=1; t<nyrs; t++) { // note this starts at 1 not 0
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

    for (int j=0; j<nyrs; j++) {
      int row_id_plus_j = row_id + j;
      out(row_id_plus_j, 0) = BoverBmsy(j);
      out(row_id_plus_j, 1) = predbio(j);
      out(row_id_plus_j, 2) = predprop(j);
      out(row_id_plus_j, 3) = predcatch(j);
      out(row_id_plus_j, 4) = ct(j);
      out(row_id_plus_j, 5) = yrs(j);
    }
    row_id = row_id + nyrs; // set up for next slot
  }
  return out;
}
