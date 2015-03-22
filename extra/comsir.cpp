#include <Rcpp.h>
using namespace Rcpp;

// comsir on first RAM stock
// priors = 60 seconds
// estimate = 120 seconds
// sumpars = 7 seconds

// [[Rcpp::export]]
DataFrame comsir_priors(int Nsim, bool logK, double K, double r, double x,
  double a, double CV, bool NormK, bool Normr, bool Norma, bool Normx,
  NumericVector Catch, double minK, double maxK, NumericVector start_r,
  bool LogisticModel, bool Obs) {

  NumericVector K_vec(Nsim);
  NumericVector r_vec(Nsim);
  NumericVector a_vec(Nsim);
  NumericVector x_vec(Nsim);
  NumericVector h(Nsim);
  NumericVector z(Nsim);
  NumericVector N1(Nsim);
  NumericVector Like(Nsim);
  NumericVector predbio(Nsim);
  NumericVector predcatch(Nsim);
  NumericVector predprop(Nsim);
  NumericVector inipredprop(Nsim);

  if (logK) {
    K_vec = exp(runif(Nsim, log(minK), log(maxK)));
  } else {
    K_vec = runif(Nsim, minK, maxK); // prior on arithmetic scale for K
  }

  if (!logK && NormK) { // normal prior on the arithmetic scale
    K_vec = rnorm(Nsim, K, CV * K);
  }

  if (Normr) {
    r_vec = rnorm(Nsim, r, CV * r);
  } else {
    r_vec = runif(Nsim, start_r(0), start_r(1)); // prior for r
  }

  if (Normx) {
    x_vec = rnorm(Nsim, x, CV * x);
  } else {
    x_vec = runif(Nsim, 0.000001, 1.0); // prior for x
  }

  if (Norma) {
    a_vec = rnorm(Nsim, a, CV * a);
  } else {
    a_vec = runif(Nsim, 0.0, 1.0); // prior for a
  }

  h = rep(0.0, Nsim); // h=0, start from unexploited state
  z = rep(1.0, Nsim); // z=1 is the Schaeffer model
  N1 = (1.0 - h) * K_vec; // initial population size

  // Like=1 if the model does not crashes, Like=0 if the model crashes
  Like = rep(1.0, Nsim);

  // Set up the numbers and project to the start of the first "real" year
  // (i.e. N1 to Nm), all simulations at once
  // initial conditions
  predbio = N1;
  predcatch = rep(Catch(0), Nsim); // the first year of catches in assumed known
  predprop = predcatch / predbio;
  inipredprop = predprop;

  for (int i=0; i<(Catch.size()-1); i++) { // TODO check this for errors

    // effort dynamics
    if (LogisticModel) { // logistic model with Bt-1
      predprop = predprop * (1.0 + x_vec * ((predbio / (a_vec * K_vec) - 1.0)));
    } else { //linear model
      predprop = predprop + x_vec * inipredprop;
    }

    // biomass dynamics
    if (Obs) {
      predbio = predbio + (r_vec * predbio * ( 1.0 - (predbio / K_vec ))) - Catch(i);
    } else {
      predbio = predbio + (r_vec * predbio * ( 1.0 - (predbio / K_vec ))) - predcatch;
    }

    predcatch = predbio * predprop;

    for (int j=0; j<Nsim; j++) {
      if ((Like(j) != 0.0 && (predbio(j) <= 0.0 || predprop(j) < 0.0 )) ||
          NumericVector::is_na(Like(j))) {
        Like(j) = 0.0;
      }
    }
  }

  // in some of the years the population crashes, we know this is impossible,
  // because the population is still there (because of the catches)

  return DataFrame::create(
      Named("N1")    = N1,
      Named("K")     = K_vec,
      Named("r")     = r_vec,
      Named("z")     = z,
      Named("a")     = a_vec,
      Named("x")     = x_vec,
      Named("h")     = h,
      Named("B")     = predbio,
      Named("Prop")  = predprop,
      Named("Like")  = Like);
}
