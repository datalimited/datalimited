#include <Rcpp.h>
using namespace Rcpp;

//' @examples
//'  ct <- c(512, 865, 511, 829, 935, 1390, 1260, 2240, 3370, 2770, 3050,
//'    3290, 4540, 3300, 3500, 3190, 2880, 3490, 5670, 6310, 9550, 8700,
//'    9130, 9160, 8490, 6400, 4420, 3680, 3190, 3960, 3290, 4220, 4220)
//'
//' out <- comsir_priors(Catch = ct,
//'    K = 800, r = 0.6, x = 0.5, a = 0.8, start_r = c(0.2, 1),
//'    minK = max(ct), maxK = max(ct) * 100, Nsim = 1e5))
// [[Rcpp::export]]
DataFrame comsir_priors(NumericVector Catch, double K, double r, double x,
  double a, NumericVector start_r,
  double minK, double maxK,
  bool logK = true,
  double CV = 0.4, bool NormK = false, bool Normr = false,
  bool Norma = false, bool Normx = false,
  bool LogisticModel = true, bool Obs = false, int Nsim = 2000L) {

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

//' @examples
//' posfun(c(0.0, 0.1, 0.2))
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

//' @examples
//'  ct <- c(512, 865, 511, 829, 935, 1390, 1260, 2240, 3370, 2770, 3050,
//'    3290, 4540, 3300, 3500, 3190, 2880, 3490, 5670, 6310, 9550, 8700,
//'    9130, 9160, 8490, 6400, 4420, 3680, 3190, 3960, 3290, 4220, 4220)
//'
//' o <- comsir_priors(Catch = ct,
//'    K = 800, r = 0.6, x = 0.5, a = 0.8, start_r = c(0.2, 1),
//'    minK = max(ct), maxK = max(ct) * 100, Nsim = 1e3)
//' o <- o[o$Like != 0, ]
//' out <- with(o, comsir_est(N1 = N1, K = K, r = r, a = a, x = x, h = h, z = z,
//'    Like = Like, ct))
// [[Rcpp::export]]
DataFrame comsir_est(NumericVector N1,
    NumericVector K,
    NumericVector r,
    NumericVector a,
    NumericVector x,
    NumericVector h,
    NumericVector z,
    NumericVector Like,
    NumericVector Catch,
    double CV = 0.4,
    bool LogisticModel = false,
    bool NormalL = true) {

  // NormallL=T #for normal likelihood
  // NormalL=F for lognormal likelihood
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

  int Nsim = K.size();
  NumericVector predbio(Nsim);
  NumericVector predcatch(Nsim);
  NumericVector predprop(Nsim);
  NumericVector inipredprop(Nsim);
  NumericVector logLike(Nsim);
  NumericVector pen1(Nsim);
  NumericVector pen2(Nsim);
  NumericMatrix temp1(Nsim, 2);
  NumericMatrix temp2(Nsim, 2);
  NumericVector CumLogLike(Nsim);
  NumericVector Mysd;
  double TotLike;

  // the parameters come from the resampling part
  // sigma = rnorm(Nsim,0.001,0.0001)//  normal informative prior for sigma
  // Set up the numbers and project to the start of the first "real" year
  // (i.e. N1 to Nm), all simulations at once
  // initial conditions
  logLike    = rep(0.0, Nsim);
  pen1       = rep(0.0, Nsim);
  pen2       = rep(0.0, Nsim);
  //temp1      = NumericMatrix(Nsim, 2);
  //temp2      = NumericMatrix(Nsim, 2);
  CumLogLike = rep(0.0, Nsim);

  for (int i=0; i<(Catch.size()-i); i++) { // TODO check this for errors
    if (i==0) {
      predbio = N1;
      predcatch = Catch(0); // first catch is assumed known -- BIG ASSUMPTION ---
      predprop = predcatch / predbio;
      inipredprop = predprop;
    } else {
      // effort dynamics
      if (LogisticModel) { // logistic model with Bt-1
        temp1 = posfun(predprop * (1 + x * ((predbio / (a * K) - 1))));
      } else { // linear model
        temp1 = posfun(predprop + x * inipredprop);
      }
      for (int j=0; j<Nsim; j++) {
        if (temp1(j, 0) > 0.99) {
          predprop(j) = 0.99;
        } else {
          predprop(j) = temp1(j, 0);
        }
        pen1(j) = pen1(j) + temp1(j, 1);
      }

      // biomass dynamics
      temp2 = posfun(predbio + (r * predbio * (1 - (predbio / K ))) - predcatch);
      for (int j=0; j<Nsim; j++) {
        predbio(j) = temp2(j, 0);
        pen2(j) = pen2(j) + temp2(j, 1);
        predcatch(j) = predbio(j) * predprop(j);
      }
    }
    //  assumption of CV=0.4 in Vasconcellos and Cochrane
    if (NormalL) {
      Mysd = CV * predcatch;
      for (int j=0; j<Nsim; j++) {
        logLike(j) = R::dnorm(Catch(i), predcatch(j), Mysd(j), 1); // normal Log likelihood
        CumLogLike(j) = CumLogLike(j) + logLike(j); // cummulative normal loglikelihood
      }
    } else {
      for (int j=0; j<Nsim; j++) {
        logLike(j) = R::dlnorm(Catch(i), log(predcatch(j)), CV, 1);
        CumLogLike(j) = CumLogLike(j) + logLike(j);
      }
    }
  }
  logLike = rep(0.0, Nsim); //  set vector to 0
  logLike = CumLogLike - pen1 - pen2; // the penalties turn out too big
  TotLike = sum(logLike);
  Like    = exp(logLike) / sum(exp(logLike)); // re-normalize so it will add up to 1

  // in some of the years the population crashes, we know this is impossible,
  // because the population is still there (because of the catches

  return DataFrame::create(
      Named("N1")       = N1,
      Named("K")        = K,
      Named("r")        = r,
      Named("z")        = z,
      Named("a")        = a,
      Named("x")        = x,
      Named("h")        = h,
      Named("B")        = predbio,    // biomass in the latest year
      Named("Prop")     = predprop,   // harvest rate in the latest year
      Named("LogLike")  = CumLogLike, // log likelihood used for DIC
      Named("Like")     = Like);      // normalized likelihood for SIR
}

// R comsir on first RAM stock
// priors = ~60 seconds
// estimate = ~120 seconds
// sumpars = ~7 seconds
