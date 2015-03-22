#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector posfun(double x, double eps) {
  // returns a new value to replace the old and a penalty for the likelihood
  NumericVector out(2);
  double y;
  double pen;
  double newvalue;

  if (x >= eps) {
    out(0) = x;
    out(1) = 0.0;
  } else {
    pen = 0.01 * pow(x - eps, 2);
    y = 1.0 - x / eps;
    newvalue = eps * (1.0 / (1.0 + y + y * y));
    out(0) = newvalue;
    out(1) = pen;
    out = c(newvalue, pen);
  }
  return out;
}

// [[Rcpp::export]]
NumericDataframe Estimate(NumericVector N1,
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

  Nsim = K.size();
  NumericVector predbio(Nsim);
  NumericVector predcatch(Nsim);
  NumericVector predprop(Nsim);
  NumericVector inipredprop(Nsim);
  NumericVector logLike(Nsim);
  NumericVector pen1(Nsim);
  NumericVector pen2(Nsim);
  NumericVector temp1(Nsim);
  NumericVector temp2(Nsim);
  NumericVector CumLogLike(Nsim);

  // the parameters come from the resampling part
  // sigma = rnorm(Nsim,0.001,0.0001)//  normal informative prior for sigma
  // Set up the numbers and project to the start of the first "real" year
  // (i.e. N1 to Nm), all simulations at once
  // initial conditions
  logLike    = rep(0.0, Nsim);
  pen1       = rep(0.0, Nsim);
  pen2       = rep(0.0, Nsim);
  temp1      = matrix(0.0, ncol = 2, nrow = Nsim); // TODO
  temp2      = matrix(0.0, ncol = 2, nrow = Nsim); // TODO
  CumLogLike = rep(0.0, Nsim);

  for (int i=0; i<(Catch.size()-i); i++) { // TODO check this for errors
    if(i==0) {
      predbio = N1;
      predcatch = Catch(0); // first catch is assumed known -- BIG ASSUMPTION ---
      predprop = predcatch / predbio;
      inipredprop = predprop;
    } else {
      // effort dynamics
      if (LogisticModel) { // logistic model with Bt-1
        temp1 = sapply(predprop * (1 + x * ((predbio / (a * K) - 1))), FUN = posfun, eps = 0.00001, simplify = TRUE);
      } else { // linear model
        temp1 = sapply((predprop + x * inipredprop), FUN = posfun, eps = 0.00001, simplify = TRUE);
        }
      predprop = pmin(temp1[1,], 0.99); // maximum harvest rate is always 0.99
      pen1 = pen1 + temp1[2, ];

      // biomass dynamics
      temp2 = sapply(predbio + (r * predbio * (1 - (predbio / K )))- predcatch,
          FUN=posfun, eps=0.00001, simplify=TRUE);
      predbio = temp2[1, ];
      pen2 = pen2 + temp2[2, ];
      predcatch = predbio * predprop;
    }
    //  assumption of CV=0.4 in Vasconcellos and Cochrane
    if (NormalL) {
      Mysd = CV * predcatch;
      logLike = dnorm(Catch(i), mean = predcatch, sd = Mysd, log = TRUE); // normal Log likelihood
      CumLogLike = CumLogLike + logLike; // cummulative normal loglikelihood
    } else {
      logLike = dlnorm(Catch(i), meanlog=log(predcatch),  sdlog=CV,  log=TRUE);
      CumLogLike = CumLogLike + logLike;
    }
  }
  logLike = rep(0.0, Nsim); //  set vector to 0
  logLike = CumLogLike - pen1 - pen2; // the penalties turn out too big
  TotLike = sum(logLike, na.rm = TRUE);
  Like    = exp(logLike) / sum(exp(logLike)); // re-normalize so it will add up to 1

  // in some of the years the population crashes, we know this is impossible,
  // because the population is still there (because of the catches

  return DataFrame::create(
      Named("N1")       = N1,
      Named("K")        = K_vec,
      Named("r")        = r_vec,
      Named("z")        = z,
      Named("a")        = a_vec,
      Named("x")        = x_vec,
      Named("h")        = h,
      Named("B")        = predbio,    // biomass in the latest year
      Named("Prop")     = predprop,   // harvest rate in the latest year
      Named("LogLike")  = CumLogLike, // log likelihood used for DIC
      Named("Like")     = Like);      // normalized likelihood for SIR
}
