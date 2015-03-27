#' COM-SIR method
#'
#' Catch-only model with sample importance resampling based on the method
#' described in Vasconcellos and Cochrane (2005).
#'
#' @param yrs A time series of years associated with the catch
#' @param ct A time series of catch
#' @param x Intrinsic rate of increase in effort. See Vasconcellos and Cochrane
#'   (2005).
#' @param k Stock size at carrying capacity
#' @param r Intrinsic rate of population growth
#' @param a Fraction of K at bioeconomic equilibrium TODO: is this correct?
#' @param start_r A numeric vector of length 2 giving the lower and upper
#'   bounds on the population growth rate parameter. This can either be
#'   specified manually or by translating resiliency categories via the function
#'   \code{\link{resilience}}.
#' @param mink Minimum possible stock size at carrying capacity
#' @param maxk Max possible stock size at carrying capacity
#' @param norm_k Logical: should \code{k} have a normal prior (\code{TRUE}) or a
#'   uniform prior between \code{mink} and \code{maxk} (\code{FALSE})?
#' @param logk If \code{logk = TRUE} and \code{norm_k = FALSE} the prior on
#'   \code{k} will be an exponentiated form a uniform distribution on a log
#'   scale between \code{log(mink)} and \code{log(maxk)}.
#' @param norm_r Logical: should {r} have a normal prior (\code{TRUE}) or a
#'   uniform prior between \code{start_r[1]} and \code{start_r[2]}
#'   (\code{FALSE})?
#' @param norm_a Logical: should \code{a} have a normal prior (\code{TRUE}) or a
#'   uniform prior between 0 and 1 (\code{FALSE})?
#' @param norm_x Logical: should \code{x} have a normal prior (\code{TRUE}) or a
#'   uniform prior between 0.000001 and 1 (\code{FALSE})?
#' @param nsim Number of iterations to run before sampling
#' @param cv Coefficient of variation on all normal prior distributions (if a
#'   given parameter is given a normal distribution instead of uniform
#'   distribution via one of the \code{norm_*} arguments.
#' @param logistic_model Logical: if \code{TRUE} then the effort dynamics model
#'   will be a logistic model. If \code{FALSE} then the effort dynamics model will
#'   be linear.
#' @param obs Logical: if \code{FALSE} then catches will be assumed to be
#'   constant from the first observation. TODO Is this right??
#' @param n_posterior Number of posterior samples to draw
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
#' @importFrom plyr adply
#' @examples
#' # TODO K and r values?
#' x <- comsir(ct = blue_gren$ct, yrs = blue_gren$yrs, k = 800, r = 0.6)
#' head(x)
#' par(mfrow = c(1, 2))
#' hist(x$BoverBmsy)
#' with(x$posterior, plot(r, k))
NULL

comsir <- function(yrs, ct, k, r, x = 0.5, a = 0.8,
  start_r = resilience(NA),
  mink = max(ct),
  maxk = max(ct) * 100, logk = TRUE, norm_k = FALSE, norm_r = FALSE,
  norm_a = FALSE, norm_x = FALSE, nsim = 2000L, cv = 0.4, logistic_model = TRUE,
  obs = FALSE, n_posterior = 1000L) {

  # TODO: note that the resilience categories were slightly different from
  # CMSY initially. Was missing 'medium'.

  o <- comsir_priors(ct = ct,
    k = k, r = r, x = x, a = a, start_r = start_r,
    mink = mink, maxk = maxk, logk = logk, norm_k = norm_k, norm_r = norm_r,
    norm_a = norm_a, norm_x = norm_x, logistic_model = logistic_model, obs = obs,
    cv = cv, nsim = nsim)
  o <- o[o$like != 0, ]

  est <- comsir_est(n1 = o$n1, k = o$k, r = o$r, a = o$a, x = o$x, h = o$h, z = o$z,
    like = o$like, ct = ct, )

  comsir_resample(est$k, est$r, est$a, est$x, est$h, est$like, yrs = yrs,
    n_posterior = n_posterior, ct, Plot = FALSE)
}

# TODO Clean up these functions!
# TODO C++ effortdyn()? Benchmark it and check
effortdyn <- function(h, k, r, x, a, yrs, ct, logistic_model) {
  # true parameters
  hrate <- h
  m <- length(ct)
  predbio<-predprop<-predcatch <- rep(0, m)

  # initial conditions
  predbio[1] <- (1.0 - hrate) * k
  predprop[1] <- ct[1]/predbio[1]
  predcatch[1] <- ct[1]

  for (t in 2:m) {
    #biomass dynamics
    predbio[t] <- predbio[t-1]+ (r*predbio[t-1] * (1-(predbio[t-1]/k))) - predcatch[t-1]
    #effort dynamics
    if (logistic_model) { # logistic model with Bt-1
      predprop[t] <- predprop[t-1]*(1+x*((predbio[t-1]/(k*a)) -1))
    }
    else { # linear model
      predprop[t] <- predprop[t-1] + x*predprop[1]
    }
    predcatch[t] <- predbio[t]*predprop[t]
  }
  BMSY <- k/2.0
  BoverBmsy <- predbio/BMSY
  xx <- cbind(BoverBmsy, predbio, predprop, predcatch, ct, yrs)
  names(xx) <- c("BoverBmsy", "biomass", "predprop", "predcatch", "obscatch", "years")
  as.data.frame(xx)
}

# @examples
# comsir_resample(k = c(100, 101, 102), h = c(0.5, 0.5, 0.5),
#   r = c(0.1, 0.2, 0.1), a = c(1, 2, 3), x = c(1, 2, 3), like = c(0, 1, 1),
#   n_posterior = 2, ct = rlnorm(10), yrs = 1:10)
comsir_resample <- function(k, r, a, x, h, like, yrs, n_posterior = 1000, ct,
  Plot = FALSE) {

  nsim <- length(k)

  # if NA the prob will be zero 20 Jan 2013
  NA.Prob <- sum(as.numeric(is.na(like)))
  ind1<-which(is.na(like))
  like[ind1]<-rep(0,length(ind1))

  # sample with replacement (just the index)
  sample_indices <- sample(nsim, n_posterior, replace = TRUE, prob = like)
  sample_indices <- sort(sample_indices)
  Post <- data.frame(h, k, r, a, x, like)[sample_indices, ]

  g <- adply(Post, 1, function(x) effortdyn(h = x$h, k = x$k, r = x$r,
    a = x$a, x = x$x, yrs = yrs, ct = ct, logistic_model = TRUE))
  g$residual <- g$ct-g$predcatch

  # diagnostics
  # TODO check if table(table()) is really what we want: (and not using this currently)
  repeated_samples <- table(table(sample_indices)) # number of sample_indices with 1, 2 or more samples
  # MSD - maximum single density,should be less than 1%
  MSD <- max(table(sample_indices)) / n_posterior * 100
  if (MSD >= 1) warning(paste0("Maximum single density was ", round(MSD, 2), "% but ",
    "should probably be < 1%."))

  MIR <- max(Post$like) / sum(Post$like) # maximum importance ratio or maximm importance weight
  cv_ir <- ((1/length(Post$like))*sum((Post$like)^2)) -
           ((1/length(Post$like))*sum(Post$like))^2
  cv_ir <- sqrt(cv_ir)/((length(Post$like)^(-0.5))*sum(Post$like))

  if (MIR >= 0.04) warning(paste0("Maximum importance ratio was ", round(MIR, 2),
      " but should probably be < 0.04."))

  if (cv_ir >= 0.04) warning(paste0("CV importance ratio was ", round(cv_ir, 2),
      " but should probably be < 0.04."))

  # new lines 06 Jan 2013  create a file with the results
  # write.table(repeated_samples, "MSD.txt")
  # write.table(MIR, "MIR.txt")
  # write.table(cv_ir, "cvIR.txt")

  # Raftery and Bao 2010 diagnostics - CHECK
  # (1) Maximum importance weight = MIR
  # (2) variance of the importance weights
  Weights <- Post$like/sum(Post$like)
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

  # TODO clean up and document the output; some can be as one data frame
  list(posterior=Post, BoverBmsy=g$BoverBmsy,Biomass=g$biomass,E=g$predprop,
    C.hat=g$predcatch,Res=g$residual,Diagno=diagnostics, MSD = MSD)
}
