#' COM-SIR method
#'
#' Catch-only model with sample importance resampling based on the method
#' described in Vasconcellos and Cochrane (2005).
#'
#' @param yr A time series of years associated with the catch
#' @param ct A time series of catch
#' @param start_r A numeric vector of length 2 giving the lower and upper
#'   bounds on the population growth rate parameter. This can either be
#'   specified manually or by translating resiliency categories via the function
#'   \code{\link{resilience}}.
#' @param k_bounds Minimum and maximum possible stock size at carrying capacity
#' @param x_bounds Minimum and maximum possible x (effort dynamics parameter) values
#' @param a_bounds Minimum and maximum possible a (effort dynamics parameter) values
#' @param norm_k Logical: should \code{k} have a normal prior (\code{TRUE}) or a
#'   uniform prior between \code{mink} and \code{maxk} (\code{FALSE})?
#' @param logk If \code{logk = TRUE} and \code{norm_k = FALSE} the prior on
#'   \code{k} will be an exponentiated form a uniform distribution on a log
#'   scale between \code{log(mink)} and \code{log(maxk)}.
#' @param nsim Number of iterations to run before sampling
#' @param logistic_model Logical: if \code{TRUE} then the effort dynamics model
#'   will be a logistic model. If \code{FALSE} then the effort dynamics model will
#'   be linear.
#' @param n_posterior Number of posterior samples to draw
#' @param obs Logical: if \code{FALSE} then a measurement-error catch model is
#'    used. If \code{TRUE} then a process-error catch model is used.
#' @param effort_bounds Lower and upper limits on rate of decrease and increase
#'   in effort from one time step to the next.
#' @param dampen Should effort dynamics parameters be excluded that may give
#'   unstable dynamics?
#' @param cv_bounds Min and max limits on residual error CV around catch
#'
#' @name comsir
#' @export
#' @references
#' Vasconcellos, M., and K. Cochrane. 2005. Overview of World Status of
#' Data-Limited Fisheries: Inferences from Landings Statistics. Pages 1-20 in G.
#' H. Kruse, V. F. Gallucci, D. E. Hay, R. I. Perry, R. M. Peterman, T. C.
#' Shirley, P. D. Spencer, B. Wilson, and D. Woodby, editors. Fisheries
#' Assessment and Management in Data-Limited Situations. Alaska Sea Grant,
#' University of Alaska Fairbanks.
#' @examples
#' # TODO K and r values?
#' x <- comsir(ct = blue_gren$ct, yr = blue_gren$yr, nsim = 1e5,
#'   n_posterior = 2e3)
#' par(mfrow = c(1, 2))
#' hist(x$quantities$bbmsy)
#' with(x$posterior, plot(r, k))
NULL

comsir <- function(yr, ct,
  start_r = resilience(NA),
  k_bounds = c(max(ct), max(ct) * 100),
  logk = TRUE,
  nsim = 1e5,
  n_posterior = 5e3, normal_like = FALSE,
  a_bounds = c(0, 1),
  x_bounds = c(1e-6, 1),
  effort_bounds = c(-Inf, Inf), obs = FALSE, dampen = FALSE,
  cv_bounds = c(0.4, 0.4)) {

  logistic_model <- TRUE
  normal_like <- FALSE

  # obs Logical: if \code{FALSE} then a measurement-error catch model is
  # used. If \code{TRUE} then a process-error catch model is used.
  # Note that this is hardcoded as FALSE because in the effort_dyn function
  # there is no obs argument

  o <- comsir_priors(ct = ct,
    start_r = start_r,
    k_bounds = k_bounds,
    logistic_model = logistic_model, obs = obs,
    nsim = nsim, a_bounds = a_bounds, x_bounds = x_bounds,
    effort_bounds = effort_bounds, dampen = dampen, cv_bounds = cv_bounds)
  o <- o[o$like > 0, ]
  # saveRDS(o, file = "prior.rds")

  est <- comsir_est(n1 = o$n1, k = o$k, r = o$r, a = o$a, x = o$x, h = o$h, z = o$z,
    like = o$like, ct = ct, logistic_model = logistic_model,
    normal_like = normal_like, cv = o$cv, obs = obs)
  # saveRDS(est, file = "est.rds")

  comsir_resample(est$k, est$r, est$a, est$x, est$h, est$like, yr = yr,
    n_posterior = n_posterior, ct, logistic_model = logistic_model, cv = est$cv)
}

# @examples
# comsir_resample(k = c(100, 101, 102), h = c(0.5, 0.5, 0.5),
#   r = c(0.1, 0.2, 0.1), a = c(1, 2, 3), x = c(1, 2, 3), like = c(0, 1, 1),
#   n_posterior = 2, ct = rlnorm(10), yr = 1:10)
comsir_resample <- function(k, r, a, x, h, like, yr, ct, n_posterior,
  logistic_model, cv) {

  nsim <- length(k)

  # if NA the prob will be zero 20 Jan 2013
  NA.Prob <- sum(as.numeric(is.na(like)))
  ind1<-which(is.na(like))
  like[ind1]<-rep(0,length(ind1))

  # sample with replacement (just the index)
  if (sum(like) > 0) {
    sample_indices <- sample(nsim, n_posterior, replace = TRUE, prob = like)
    post <- data.frame(h, k, r, a, x, like, cv)[sample_indices, ]
  } else {
    warning("All likelihoods were too small. Returning NULL.")
    return(NULL)
  }

  quantities <- effortdyn(h = post$h, k = post$k, r = post$r, a = post$a,
    x = post$x, yr = yr, ct = ct, logistic_model = logistic_model)
  quantities <- setNames(as.data.frame(quantities),
    nm = c("bbmsy", "predbio", "predprop", "predcatch", "ct", "yr"))
  quantities$residual <- quantities$ct - quantities$predcatch
  quantities$sample_id <- rep(seq_along(sample_indices), each = length(yr))

  # diagnostics
  # check if table(table()) is really what we want: (and not using this currently)
  # repeated_samples <- table(table(sample_indices)) # number of sample_indices with 1, 2 or more samples
  MSD <- max(table(sample_indices)) / n_posterior * 100
  MIR <- max(post$like) / sum(post$like) # maximum importance ratio or maximm importance weight
  cv_ir <- ((1/length(post$like))*sum((post$like)^2)) -
    ((1/length(post$like))*sum(post$like))^2
  cv_ir <- sqrt(cv_ir)/((length(post$like)^(-0.5))*sum(post$like))

  if (MSD >= 1) warning(paste0("Maximum single density was ", round(MSD, 2), "% but ",
    "should probably be < 1%."))
  if (MIR >= 0.04) warning(paste0("Maximum importance ratio was ", round(MIR, 2),
    " but should probably be < 0.04."))
  tryCatch({
    if (cv_ir >= 0.04) warning(paste0("CV importance ratio was ", round(cv_ir, 2),
      " but should probably be < 0.04."))
  }, error = function(e) warning(e))

  # Raftery and Bao 2010 diagnostics - CHECK
  # (1) Maximum importance weight = MIR
  # (2) variance of the importance weights
  Weights <- post$like/sum(post$like)
  Var.RW <- (n_posterior * Weights -1)^2 #vector
  Var.RW <- (1 / n_posterior) * sum(Var.RW) #scalar
  # (3) entropy of the importance weights relative to uniformity
  Entropy <- -1 * sum(Weights * (log(Weights) / log(n_posterior)))
  # (4)Expected number of unique points after resampling
  Exp.N <- sum((1 - (1 - Weights)^n_posterior))
  # Effective sample size ESS
  ESS <- 1 / (sum(Weights^2))

  diagnostics <- c(nsim, NA.Prob, MIR, cv_ir, MSD, Var.RW, Entropy, Exp.N, ESS)
  names(diagnostics) <- c("N.nocrash.priors", "NA.Prob", "MIR", "cv_ir", "MSD",
    "Var.RW", "Entropy", "Exp.N", "ESS")

  bbmsy <- reshape2::dcast(quantities, sample_id ~ yr, value.var = "bbmsy")[, -1]
  bbmsy <- data.frame(year = yr, catch = ct, summarize_bbmsy(bbmsy))

  list(posterior = post, quantities = quantities, diagnostics = diagnostics,
    msd = MSD, bbmsy = bbmsy)
}

comsir_priors <- function(ct, start_r,
  nsim, x_bounds, a_bounds, k_bounds, effort_bounds, logistic_model, obs, dampen,
  cv_bounds) {

  priors <- comsir_bounds(x_bounds = x_bounds, a_bounds = a_bounds,
    k_bounds = k_bounds, effort_bounds = effort_bounds, nsim = nsim, dampen = dampen)
  colnames(priors) <- c("a", "x", "k", "effort_bt0", "effort_btk")

  k_vec <- priors[,"k"]
  a_vec <- priors[,"a"]
  x_vec <- priors[,"x"]
  r_vec <- runif(nsim, start_r[1], start_r[2])
  cv_vec <- runif(nsim, cv_bounds[1], cv_bounds[2])

  h <- rep(0, nsim)
  z <- rep(1, nsim)
  n1 <- (1 - h) * k_vec

  like <- rep(1, nsim)
  predbio <- n1
  predcatch <- ct[1]
  predprop <- predcatch / predbio
  inipredprop <- predprop
  m <- length(ct)

  # for debugging (save catch series)
  # check <- matrix(nrow = m, ncol = length(k_vec))

  for (i in 1:m) {
    if (logistic_model)
      predprop <- predprop * (1 + x_vec * ((predbio/(a_vec * k_vec) - 1)))
    else predprop = predprop + x_vec * inipredprop
    if (obs)
      predbio = predbio + (r_vec * predbio * (1 - (predbio/k_vec))) - ct[i]
    else predbio = predbio + (r_vec * predbio * (1 - (predbio/k_vec))) - predcatch
    predcatch = predbio * predprop

    # check[i, ] <- predcatch
    # can this just be vectorized?
    like <- check_comsir_lik(nsim, predbio, predprop, like)
    #for (i in 1:nsim) {
    #if ((like[i] != 0 & (predbio[i] <= 0 || predprop[i] < 0)) ||
    #is.na(like[i])) like[i] <- 0
    #}
  }

  data.frame(n1, k = k_vec, r = r_vec, z, a = a_vec, x = x_vec, h, biomass =
      predbio, prop = predprop, like = like, cv = cv_vec)
}

comsir_est <- function(n1, k, r, a, x, h, z, like, ct, cv,
  logistic_model, normal_like, obs) {

  #  normal_like=true for normal likelihood
  #  normal_like=false for lognormal likelihood
  #
  #  The log normal distribution has density
  #   f(x) = 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))
  #      where mu and sigma are the mean and standard deviation of the
  #      logarithm. The mean is E(X) = exp(mu + 1/2 sigma^2), and the
  #      variance Var(X) = exp(2*mu + sigma^2)*(exp(sigma^2) - 1) and hence
  #      the coefficient of variation is sqrt(exp(sigma^2) - 1) which is
  #      approximately sigma when that is small (e.g., sigma < 1/2).
  #     x<-seq(0.001,100,0.001)
  #     y<-dlnorm(log(x),meanlog=log(1),sdlog=1,log=FALSE)
  #      dlnorm(1) == dnorm(0)
  #     z<-dnorm(log(x),mean=0,sdlog=1,log=FALSE)
  #     plot(log(x),y,type="l")
  #     plot(log(x),z,type="l")

  nsim = length(k) #.size();

  # the parameters come from the resampling part
  # sigma = rnorm(nsim,0.001,0.0001)//  normal informative prior for sigma
  # Set up the numbers and project to the start of the first "real" year
  # (i.e. n1 to Nm), all simulations at once
  # initial conditions
  loglike     = rep(0.0, nsim)
  pen1        = rep(0.0, nsim)
  pen2        = rep(0.0, nsim)
  cum_loglike = rep(0.0, nsim)

  for (i in seq_along(ct)) {
    if (i==1) {
      predbio = n1
      predcatch = rep(ct[1], nsim) # first catch is assumed known
      predprop = predcatch / predbio
      inipredprop = predprop
    } else {
      # effort dynamics
      if (logistic_model) { # logistic model with Bt-1
        temp1 = posfun(predprop * (1 + x * ((predbio / (a * k) - 1))))
      } else { # linear model
        temp1 = posfun(predprop + x * inipredprop)
      }

      predprop[temp1[,1] >  0.99] <- 0.99
      predprop[temp1[,1] <= 0.99] <- temp1[temp1[,1] <= 0.99, 1]
      pen1 = pen1 + temp1[, 2]

      # biomass dynamics
      if (obs) {
        temp2 = posfun(predbio + (r * predbio * (1 - (predbio/k))) - ct[i])
      } else {
        temp2 = posfun(predbio + (r * predbio * (1 - (predbio / k ))) - predcatch)
      }
      predbio = temp2[, 1]
      pen2 = pen2 + temp2[, 2]
      predcatch = predbio * predprop
    }
    #  assumption of cv=0.4 in Vasconcellos and Cochrane
    if (normal_like) {
      my_sd = cv * predcatch
      loglike = dnorm(ct[i], predcatch, my_sd, log = TRUE) # normal Log likelihood
      cum_loglike = cum_loglike + loglike # cummulative normal loglikelihood
    } else {
      loglike = dlnorm(ct[i], log(predcatch), cv, log = TRUE)
      cum_loglike = cum_loglike + loglike;
    }
    # print(predbio[1])
    # print(c(predbio,predprop,predcatch,ct[i]))
  }
  loglike = rep(0.0, nsim) #  set vector to 0
  loglike = cum_loglike - pen1 - pen2 # the penalties turn out too big
  Totlike = sum(loglike)
  # if (length(unique(exp(loglike))))
  # like    = exp(loglike) / sum(exp(loglike)) # re-normalize so it will add up to 1
  # There can be underflow issues here
  # One possible solution:
  # http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability/66621#66621
  # Another solution (but very slow)
  # loglike2 <- Rmpfr::mpfr(loglike, prec = 25L) # prevent underflow
  # like    = as.numeric(exp(loglike2) / sum(exp(loglike2)))
  like = as.numeric(exp(loglike) / sum(exp(loglike))) # re-normalize so it will add up to 1

  # in some of the years the population crashes, we know this is impossible,
  # because the population is still there (because of the catches

  data.frame(n1, k, r, z, a, x, h, cv,
    B        = predbio,    # biomass in the latest year
    prop     = predprop,   # harvest rate in the latest year
    loglike  = cum_loglike, # log likelihood used for DIC
    like     = like)      # normalized likelihood for SIR
}
