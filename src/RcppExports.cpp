// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// schaefer_cmsy
DataFrame schaefer_cmsy(NumericVector r_lim, NumericVector k_lim, double sig_r, NumericVector startbio, NumericVector yr, NumericVector ct, int interyr_index, double prior_log_mean, double prior_log_sd, NumericVector interbio, int reps);
RcppExport SEXP datalimited_schaefer_cmsy(SEXP r_limSEXP, SEXP k_limSEXP, SEXP sig_rSEXP, SEXP startbioSEXP, SEXP yrSEXP, SEXP ctSEXP, SEXP interyr_indexSEXP, SEXP prior_log_meanSEXP, SEXP prior_log_sdSEXP, SEXP interbioSEXP, SEXP repsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type r_lim(r_limSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type k_lim(k_limSEXP );
        Rcpp::traits::input_parameter< double >::type sig_r(sig_rSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type startbio(startbioSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type yr(yrSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ct(ctSEXP );
        Rcpp::traits::input_parameter< int >::type interyr_index(interyr_indexSEXP );
        Rcpp::traits::input_parameter< double >::type prior_log_mean(prior_log_meanSEXP );
        Rcpp::traits::input_parameter< double >::type prior_log_sd(prior_log_sdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type interbio(interbioSEXP );
        Rcpp::traits::input_parameter< int >::type reps(repsSEXP );
        DataFrame __result = schaefer_cmsy(r_lim, k_lim, sig_r, startbio, yr, ct, interyr_index, prior_log_mean, prior_log_sd, interbio, reps);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// comsir_priors
DataFrame comsir_priors(NumericVector Catch, double K, double r, double x, double a, NumericVector start_r, double minK, double maxK, bool logK = true, double CV = 0.4, bool NormK = false, bool Normr = false, bool Norma = false, bool Normx = false, bool LogisticModel = true, bool Obs = false, int Nsim = 2000L);
RcppExport SEXP datalimited_comsir_priors(SEXP CatchSEXP, SEXP KSEXP, SEXP rSEXP, SEXP xSEXP, SEXP aSEXP, SEXP start_rSEXP, SEXP minKSEXP, SEXP maxKSEXP, SEXP logKSEXP, SEXP CVSEXP, SEXP NormKSEXP, SEXP NormrSEXP, SEXP NormaSEXP, SEXP NormxSEXP, SEXP LogisticModelSEXP, SEXP ObsSEXP, SEXP NsimSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type Catch(CatchSEXP );
        Rcpp::traits::input_parameter< double >::type K(KSEXP );
        Rcpp::traits::input_parameter< double >::type r(rSEXP );
        Rcpp::traits::input_parameter< double >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type a(aSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type start_r(start_rSEXP );
        Rcpp::traits::input_parameter< double >::type minK(minKSEXP );
        Rcpp::traits::input_parameter< double >::type maxK(maxKSEXP );
        Rcpp::traits::input_parameter< bool >::type logK(logKSEXP );
        Rcpp::traits::input_parameter< double >::type CV(CVSEXP );
        Rcpp::traits::input_parameter< bool >::type NormK(NormKSEXP );
        Rcpp::traits::input_parameter< bool >::type Normr(NormrSEXP );
        Rcpp::traits::input_parameter< bool >::type Norma(NormaSEXP );
        Rcpp::traits::input_parameter< bool >::type Normx(NormxSEXP );
        Rcpp::traits::input_parameter< bool >::type LogisticModel(LogisticModelSEXP );
        Rcpp::traits::input_parameter< bool >::type Obs(ObsSEXP );
        Rcpp::traits::input_parameter< int >::type Nsim(NsimSEXP );
        DataFrame __result = comsir_priors(Catch, K, r, x, a, start_r, minK, maxK, logK, CV, NormK, Normr, Norma, Normx, LogisticModel, Obs, Nsim);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// posfun
NumericMatrix posfun(NumericVector x, double eps = 0.00001);
RcppExport SEXP datalimited_posfun(SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type eps(epsSEXP );
        NumericMatrix __result = posfun(x, eps);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// comsir_est
DataFrame comsir_est(NumericVector N1, NumericVector K, NumericVector r, NumericVector a, NumericVector x, NumericVector h, NumericVector z, NumericVector Like, NumericVector Catch, double CV = 0.4, bool LogisticModel = false, bool NormalL = true);
RcppExport SEXP datalimited_comsir_est(SEXP N1SEXP, SEXP KSEXP, SEXP rSEXP, SEXP aSEXP, SEXP xSEXP, SEXP hSEXP, SEXP zSEXP, SEXP LikeSEXP, SEXP CatchSEXP, SEXP CVSEXP, SEXP LogisticModelSEXP, SEXP NormalLSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type N1(N1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type K(KSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Like(LikeSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type Catch(CatchSEXP );
        Rcpp::traits::input_parameter< double >::type CV(CVSEXP );
        Rcpp::traits::input_parameter< bool >::type LogisticModel(LogisticModelSEXP );
        Rcpp::traits::input_parameter< bool >::type NormalL(NormalLSEXP );
        DataFrame __result = comsir_est(N1, K, r, a, x, h, z, Like, Catch, CV, LogisticModel, NormalL);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
