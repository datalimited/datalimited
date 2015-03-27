effortdyn_r <- function(h, k, r, x, a, yrs, ct, logistic_model) {
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

library(dplyr)
set.seed(1)
ct <- rlnorm(10)
o1_lF <- effortdyn_r(h = 0.5, k = 100, r = 0.4, x = 0.5, a = 0.8, yrs = 1:10, ct = ct, logistic_model = FALSE)
o1_lT <- effortdyn_r(h = 0.5, k = 100, r = 0.4, x = 0.5, a = 0.8, yrs = 1:10, ct = ct, logistic_model = TRUE)

# C++ version:
o1_lF_cpp <- effortdyn(h = 0.5, k = 100, r = 0.4, x = 0.5, a = 0.8, yrs = 1:10, ct = ct, logistic_model = FALSE) %>% as.data.frame
names(o1_lF_cpp) <- names(o1_lF)

o1_lT_cpp <- effortdyn(h = 0.5, k = 100, r = 0.4, x = 0.5, a = 0.8, yrs = 1:10, ct = ct, logistic_model = TRUE) %>% as.data.frame
names(o1_lT_cpp) <- names(o1_lT)

identical(o1_lF, o1_lF_cpp)
identical(o1_lT, o1_lT_cpp)

N <- 1e2
r_fun <- function() plyr::rdply(N, function(x) effortdyn_r(h = 0.5, k = 100, r = 0.4, x = 0.5, a = 0.8, yrs = 1:10, ct = ct, logistic_model = FALSE))

cpp_fun <- function() effortdyn(h = rep(0.5, N), k = rep(100, N), r = rep(0.4, N), x = rep(0.5, N), a = rep(0.8, N), yrs = 1:10, ct = ct, logistic_model = TRUE)

library(microbenchmark)
mb <- microbenchmark(r = r_fun(), cpp = cpp_fun(), times = 100L)
mb
# about 50 times faster
