# TODO Clean up these functions!
effortdyn <- function(h, K, r, x, a, yrs, Catch, LogisticModel) {
  # true parameters
  hrate <- h
  m <- length(Catch)
  predbio<-predprop<-predcatch <- rep(0, m)

  # initial conditions
  predbio[1] <- (1.0 - hrate) * K
  predprop[1] <- Catch[1]/predbio[1]
  predcatch[1] <- Catch[1]

  for (t in 2:m) {
    #biomass dynamics
    predbio[t] <- predbio[t-1]+ (r*predbio[t-1] * (1-(predbio[t-1]/K))) - predcatch[t-1]
    #effort dynamics
    if (LogisticModel) { # logistic model with Bt-1
      predprop[t] <- predprop[t-1]*(1+x*((predbio[t-1]/(K*a)) -1))
    }
    else { # linear model
      predprop[t] <- predprop[t-1] + x*predprop[1]
    }
    predcatch[t] <- predbio[t]*predprop[t]
  }
  BMSY <- K/2.0
  BoverBmsy <- predbio/BMSY
  xx <- cbind(BoverBmsy, predbio, predprop, predcatch, Catch, yrs)
  names(xx) <- c("BoverBmsy", "biomass", "predprop", "predcatch", "obscatch", "years")
  as.data.frame(xx)
}

#' @examples
#' comsir_resample(K = c(100, 101, 102), h = c(0.5, 0.5, 0.5),
#'   r = c(0.1, 0.2, 0.1), a = c(1, 2, 3), x = c(1, 2, 3), Like = c(0, 1, 1),
#'   Npost = 2, Catch = rlnorm(10), yrs = 1:10)
comsir_resample<- function(K, r, a, x, h, Like, yrs, Npost = 1000, Catch, Plot = FALSE) {
  Nsim <- length(K)

  # if NA the prob will be zero 20 Jan 2013
  # print("how many probabilities are NAs?")
  # print(table(is.na(ParVals$Like))) # debug
  NA.Prob <- sum(as.numeric(is.na(Like)))
  ind1<-which(is.na(Like))
  Like[ind1]<-rep(0,length(ind1))

  # sample with replacement (just the index)
  Ind <- sample(Nsim, Npost, replace = TRUE, prob = Like)
  Ind <- sort(Ind)
  Post <- data.frame(h, K, r, a, x, Like)[Ind, ]

  g <- plyr::adply(Post, 1, function(x) effortdyn(h = x$h, K = x$K, r = x$r,
    a = x$a, x = x$x, yrs = yrs, Catch = Catch, LogisticModel = TRUE))
  g$residual <- g$Catch-g$predcatch

  #diagnostics
  # TODO check if table(table()) is really what we want:
  repeatedSamples <- table(table(Ind)) #number of Ind with 1, 2 or more samples
  #MSD - maximum single density,should be less than 1%
  MSD <- max(table(Ind)) / Npost * 100
  if (MSD >= 1) warning(paste0("Maximum single density was ", MSD, "% but ",
    "should probably be < 1%."))

  ##novo
  MIR <- max(Post$Like) / sum(Post$Like) ###maximum importance ratio or maximm importance weight
  CV.IR <- ((1/length(Post$Like))*sum((Post$Like)^2)) - ((1/length(Post$Like))*sum(Post$Like))^2
  CV.IR <- sqrt(CV.IR)/((length(Post$Like)^(-0.5))*sum(Post$Like))
  print(repeatedSamples)

  if (MIR >= 0.04) warning(paste0("Maximum importance ratio was ", MIR, " but ",
    "should probably be < 0.04."))

  #print("CV.IR - CV importance ratio, should be less than 0.04")
  if (CV.IR >= 0.04) warning(paste0("CV importance ratio was ", CV.IR, " but ",
    "should probably be < 0.04."))


  ##new lines 06 Jan 2013  create a file with the results
  #write.table(repeatedSamples, "MSD.txt")
  #write.table(MIR, "MIR.txt")
  #write.table(CV.IR, "CVIR.txt")
  ###end
  #Raftery and Bao 2010 diagnostics - CHECK
  #(1) Maximum importance weight = MIR
  # (2) variance of the importance weights
  Weights <- Post$Like/sum(Post$Like)
  Var.RW <- (Npost*Weights -1)^2 #vector
  Var.RW <- (1/Npost)*sum(Var.RW) #scalar
  # (3) entropy of the importance weights relative to uniformity
  Entropy <- -1* sum(Weights*(log(Weights)/log(Npost)))
  # (4)Expected number of unique points after resampling
  Exp.N <- sum((1-(1-Weights)^Npost))
  # Effective sample size ESS
  ESS <- 1/(sum(Weights^2))

  diagnostics <- c(Nsim, NA.Prob, MIR, CV.IR, MSD, Var.RW, Entropy, Exp.N, ESS)
  names(diagnostics) <- c("N.nocrash.priors", "NA.Prob", "MIR", "CV.IR", "MSD",
    "Var.RW", "Entropy", "Exp.N", "ESS")

  # TODO return the data frame part as one big data frame:
  list(posterior=Post, BoverBmsy=g$BoverBmsy,Biomass=g$biomass,E=g$predprop,C.hat=g$predcatch,Res=g$residual,Diagno=diagnostics, MSD = MSD)
}
