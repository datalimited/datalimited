runThor2<-function(){
  library(R2jags)
  library(matlab)
  RunFile = "./"
  res.list<-list()
  ##StockName = names(input)[[1]]
  stock_id<-unique(input$stockid)
  seed<-stockNumber+1
  set.seed(seed)
  ##
  stock.dat<-subset(input, stockid==stock_id[stockNumber] & !is.na(CtoUse))
  stock.dat<-subset(stock.dat, FAO.Fishing.area==unique(FAO.Fishing.area)[1]) ## only one FAO fishing area
  stock.dat<-stock.dat[order(stock.dat$tsyear),]
  Year   <- stock.dat$tsyear
  Ct   <- stock.dat$CtoUse  ## assumes that catch is given in tonnes, transforms to '000 tonnes
  StockName = stock_id[stockNumber]
  ##Ct = input[[StockName]][[1]]$catch$data
  ##Ct = Load[,'catch.data']
  MaxCt = max(Ct)
  Ct = Ct / MaxCt
  ##Year = input[[StockName]][[1]]$catch$year
  ##Year = Load[,'catch.year']
  Nyears = length( Year )
  Bt = rep(1,Nyears)
  B0 = Bt[1]
  DepletionPriorSD = 1000
  ##
  ## R prior
  ## if(input[[StockName]][[1]]$tmat<1) rPrior1 = 1
  ## if(input[[StockName]][[1]]$tmat>=1 & input[[StockName]][[1]]$tmat<4) rPrior1 = 2
  ## if(input[[StockName]][[1]]$tmat>=4 & input[[StockName]][[1]]$tmat<10) rPrior1 = 3
  ## if(input[[StockName]][[1]]$tmat>=10) rPrior1 = 4
  ## if(input[[StockName]][[1]]$tmax<3) rPrior2 = 1
  ## if(input[[StockName]][[1]]$tmax>=3 & input[[StockName]][[1]]$tmax<10) rPrior2 = 2
  ## if(input[[StockName]][[1]]$tmax>=10 & input[[StockName]][[1]]$tmax<30) rPrior2 = 3
  ## if(input[[StockName]][[1]]$tmax>=30) rPrior2 = 4
  ## ##rPrior1 = 2; rPrior2 = 2
  ## rPrior = max(rPrior1, rPrior2)
  ## if(rPrior==1){ r_min=0.6; r_max=1.5 }
  ## if(rPrior==2){ r_min=0.2; r_max=1.0 }
  ## if(rPrior==3){ r_min=0.05; r_max=0.5 }
  ## if(rPrior==4){ r_min=0.015; r_max=0.1 }
  ## CM: altered Dec 2014
  res  <- unique(stock.dat$res) ## resilience from FishBase
  start_r  <- if(res == "Very low"){c(0.015, 0.1)} else if(res == "Low") {c(0.05,0.5)} else if(res == "High") {c(0.6,1.5)}  else {c(0.2,1)} ## Medium, or default if no re
  r_min=start_r[1]
  r_max=start_r[2]
# JAGS model 1 -- Latent Bt, Process errors for Pt and Et
JagsModel = "
model {
  SigmaE ~ dunif(0.01,1)
  SigmaP <- SigmaE
  SigmaQ <- SigmaP

  #logitB1overB0 ~ dunif(-4.6,4.6)
  #lnr ~ dnorm(rPrior[1],rPrior[2])
  lnr <- log(r)
  lnB0 ~ dunif(-4.6,4.6)
  x ~ dunif(0.01,0.5)
  a ~ dunif(0.1,2)

  #r <- exp(lnr)
  #r ~ dunif(r_min,r_max)
  #MSY <- B0 * r / 4
  B0 <- exp(lnB0)
  ln_Final_Depletion <- log(Bt[Nyears]) - lnB0
  Final_Depletion <- exp(ln_Final_Depletion)
  Final_Depletion_Error_Ratio <- exp(ln_Final_Depletion) / DepletionTrue[Nyears]
  #Final_Depletion_Prior[1] ~ dlnorm( ln_Final_Depletion, Final_Depletion_Prior[2])
  MSY_min <- B0 * r_min / 4
  MSY_max <- B0 * r_max / 4
  MSY ~ dunif(MSY_min,MSY_max)
  r <- 4 * MSY / B0

  #Bt[1] <- B0 * 1/(1+exp(-logitB1overB0))
  Bt_rel[1] <- 1
  Bt[1] <- Bt_rel[1] * B0
  lnE0 ~ dunif(-11.5,0)    # ln(0.00001) - ln(0.1)
  E0 <- exp(lnE0)
  Et[1] <- E0
  Dummy <- lnE0
  Bupper <- 2*B0   # Necessary for when r is low, so that Bt doesn't drift to far above B0
  Ct_hat[1] <- 0
  for(YearI in 2:Nyears){
    # Define time-varying precision
    TauE[YearI] <- pow(SigmaE,-2) * pow(EffortSD[YearI],-2)
    TauB[YearI] <- pow(SigmaP,-2) * pow(BioSD[YearI],-2)
    TauQ[YearI] <- pow(SigmaQ,-2) * pow(Q_SD[YearI],-2)
    # Stochastic draw for Bt given Bt_exp
    Pt[YearI-1] <- r*Bt_rel[YearI-1]*B0 * ( 1 - Bt_rel[YearI-1] )
    #ln_Bt_rel_exp[YearI] <- log( max( Bt_rel[YearI-1] + Pt[YearI-1] - Ct[YearI-1], -0.0001/(Bt_rel[YearI-1] + Pt[YearI-1] - Ct[YearI-1]) ) )  # 1e-10 is about the lowest I can set this
    ln_Bt_rel_exp[YearI] <- log( max( (Bt_rel[YearI-1]*B0 + Pt[YearI-1] - Ct[YearI-1])/B0, 1e-12 ) )  # 1e-10 is about the lowest I can set this
    Bt_rel[YearI] ~ dlnorm( ln_Bt_rel_exp[YearI], TauB[YearI]) T(0.001,Bupper)
    Bt[YearI] <- Bt_rel[YearI] * B0
    # Set up next effort computation used in next year's stochastic draw for Ct
    Et[YearI] <- min( Ct[YearI] / (Bt_rel[YearI]*B0), 0.99 )
    # Stochastic draw for Ct given Ct_hat
    Ct_hat[YearI] <- Et[YearI-1] * ( Bt_rel[YearI-1] / (a/2) )^x * Bt_rel[YearI]*B0
    ln_Ct_hat[YearI] <- log( Ct_hat[YearI] )
    Ct[YearI] ~ dlnorm(ln_Ct_hat[YearI],TauQ[YearI])
  }
}
"
  ## Run Settings
  ## which effort dynamics model
  ## ModelType = 1 : Process errors in biomass and effort dynamics
  ## ModelType = 2 : same plus variability in catch equation
  ##ModelType = 1 # 1=Bt; 2=Bt+Et
  ## JAGS settings
  NburninPrelim = 1e3
  NiterPrelim = 2e3
  NthinPrelim = 1e0
  NchainsPrelim = 100
  NburninJags = 1e6
  NiterJags = 3e6
  NthinJags = 1e3
  ##NburninJags = 1e3
  ##NiterJags = 3e3
  ##NthinJags = 1e1
  Nchains = 3
  ## Note! ModelType 1
  ModelType <- 1
  ## Write JAGS model
  cat(JagsModel, file=paste(RunFile,"dnorm.bug",sep=""))
  Params2Save = c("Ct_hat","Pt","Bt","B0","r","a","x","SigmaE","SigmaP","MSY","Final_Depletion","Final_Depletion_Error_Ratio","E0") #,"logitB1overB0")
  ## Run JAGS
  ## Nyears = number of years
  ## Ct = catch time series
  ## Q_SD = relative index of CV in catch equation (i.e. catchability)
  ## BioSD = relative index of CV in biomass process errors
  ## EffortSD = relative index of CV in effort dynamics process errors
  ## DepletionTrue = true Bt/B0 (only used for computing errors)
  ## Final_Depletion_Prior = mean and precision for prior on final Bt/B0
  ## r_min and r_max = min and max for uniform prior on r (intrinsic growth rate)
  Seed = ceiling(runif(1,0,1e6))
  ##Seed <- 111006
  set.seed(Seed)
  ## Final_Depletion_Prior=c(Bt[Nyears]/B0,1/DepletionPriorSD^2),
  DataJags = list(Nyears=Nyears, Ct=Ct, Q_SD=rep(1,Nyears), BioSD=rep(1,Nyears), EffortSD=rep(1,Nyears), DepletionTrue=Bt/B0, r_min=r_min, r_max=r_max) #, SigmaE=SigmaE, SigmaP=SigmaP)
  ##
  inits_prelim = function(){
    lnE0 = log(runif(1, 0.00002, 0.999999))
    x = runif(1, 0.01,0.5)
    a = runif(1, 0.1,2)
    r = runif(1, r_min,r_max)
    lnB0 = log(max(Ct)/r*4) + runif(1, log(0.2), log(25))
    MSY = r * exp(lnB0) / 4
    Bt_rel = runif( length(Ct), min=Ct/exp(lnB0), max=1 )
    Bt_rel[1] <- NA
    SigmaE = runif(1, 0.01,1)
    Return = list( lnE0=lnE0, x=x, a=a, MSY=MSY, lnB0=lnB0, Bt_rel=Bt_rel, SigmaE=SigmaE )
    ##Return = list( lnE0=lnE0, x=x, a=a, r=r, lnB0=lnB0, Bt_rel=Bt_rel, SigmaE=SigmaE )
    return(Return)
  }
  inits_resample = function(){
    Which = sample(size=1, 1:length(RelLike), prob=RelLike)
    lnE0 = log(Jags_prelim$BUGSoutput$sims.list$E0)[Which]
    x = Jags_prelim$BUGSoutput$sims.list$x[Which]
    a = Jags_prelim$BUGSoutput$sims.list$a[Which]
    r = Jags_prelim$BUGSoutput$sims.list$r[Which]
    MSY = Jags_prelim$BUGSoutput$sims.list$MSY[Which]
    lnB0 = log(Jags_prelim$BUGSoutput$sims.list$B0)[Which]
    Bt_rel = Jags_prelim$BUGSoutput$sims.list$Bt[Which,] / Jags_prelim$BUGSoutput$sims.list$B0[Which]
    Bt_rel[1] <- NA
    SigmaE = Jags_prelim$BUGSoutput$sims.list$SigmaE[Which]
    Return = list( lnE0=lnE0, x=x, a=a, MSY=MSY, lnB0=lnB0, Bt_rel=Bt_rel, SigmaE=SigmaE )
    ##Return = list( lnE0=lnE0, x=x, a=a, r=r, lnB0=lnB0, Bt_rel=Bt_rel, SigmaE=SigmaE )
    return(Return)
  }
  ##
  ##jags.model(file=paste(RunFile,"dnorm.bug",sep=""), data=DataJags, inits=inits_prelim, n.chains=NchainsPrelim, quiet=T)
  ##update(dm.mod, n.iter=mcmc.burn)
  ##Samples = coda.samples(dm.mod, variable.names=pars.to.save, n.iter=mcmc.chainLength, thin=mcmc.thin)
  ## Try to initialize the model
  Iteration = IterationPrelim = 1
  tic()
  while(IterationPrelim < 1000){
    ##
    NchainsIteration = ceiling(NchainsPrelim * (1000-2*IterationPrelim)/1000)   # Spend half of attempts at Nchains = 3
    if(NchainsIteration<1) NchainsIteration = 1
    Try = try(Jags_prelim <- jags(model.file=paste(RunFile,"dnorm.bug",sep=""), inits=inits_prelim, working.directory=NULL, data=DataJags, parameters.to.save=Params2Save, n.chains=NchainsIteration, n.thin=NthinPrelim, n.iter=NiterPrelim, n.burnin=NburninPrelim))
    IterationPrelim = ifelse(class(Try)!='try-error', 1e20, IterationPrelim+1)
    print(IterationPrelim)
  }
  ## If it initializes, then try to run the model
  if(IterationPrelim==1e20){
    RelLike = exp(-(Jags_prelim$BUGSoutput$sims.list$deviance-min(Jags_prelim$BUGSoutput$sims.list$deviance)))
    ##
    while(Iteration < 1000){
      Try = try( Jags <- jags(model.file=paste(RunFile,"dnorm.bug",sep=""), inits=inits_resample, working.directory=NULL, data=DataJags, parameters.to.save=Params2Save, n.chains=Nchains, n.thin=NthinJags, n.iter=NiterJags, n.burnin=NburninJags) )
      Iteration = ifelse(class(Try)!='try-error', 1e20, Iteration+1)
    }
  }
  Time = toc(echo=FALSE)
  ##
  ## If it both initializes and runs, then save results
  if(IterationPrelim==1e20 & Iteration==1e20){
    BoverBmsy = Jags$BUGSoutput$sims.list$Bt / outer(as.vector(Jags$BUGSoutput$sims.list$B0),rep(1,Nyears)) * 2
    Neff = Jags$BUGSoutput$summary[,'n.eff']
    Neff = ifelse(Neff==1,NA,Neff)
    Neff = min(Neff, na.rm=TRUE)
    if(ModelType==1) E = array(0, dim=dim(Jags$BUGSoutput$sims.list$Bt))
    if(ModelType==2) E = Jags$BUGSoutput$sims.list$Et
  }else{
    BoverBmsy = matrix(0, ncol=Nyears, nrow=(NiterJags-NburninJags)/NthinJags*Nchains)
    Neff = rep(3, 100)
    names(Neff)[1:2] = c("Ct_hat[1]","Pt[1]")
    Neff = Neff[-match(c("Ct_hat[1]","Pt[1]"),names(Neff))]
    Neff = min(Neff)
    E = array(0, dim=dim(BoverBmsy))
  }
  Results = suppressWarnings(data.frame('stock_id'=StockName, 'b_bmsy'=apply(BoverBmsy,MARGIN=2,FUN=mean), 'b_bmsyUpper'=apply(BoverBmsy,MARGIN=2,FUN=quantile,prob=0.925), 'b_bmsyLower'=apply(BoverBmsy,MARGIN=2,FUN=quantile,prob=0.075), 'b_bmsy_iq25'=apply(BoverBmsy,MARGIN=2,FUN=quantile,prob=0.25), 'b_bmsy_iq75'=apply(BoverBmsy,MARGIN=2,FUN=quantile,prob=0.75), 'E'=apply(E,MARGIN=2,FUN=mean), 'E_Upper'=apply(E,MARGIN=2,FUN=quantile,prob=0.925), 'E_Lower'=apply(E,MARGIN=2,FUN=quantile,prob=0.075), 'E_iq25'=apply(E,MARGIN=2,FUN=quantile,prob=0.25), 'E_iq75'=apply(E,MARGIN=2,FUN=quantile,prob=0.75), 'year'=Year, 'seed'=Seed, 'convergence'=ifelse(Neff>200,"Strong",ifelse(Neff>30,"Weak","Not")), 'n_iterations'=NiterJags*Nchains, 'effective_sample_size'=Neff, 'run_time'=Time['elapsed'], 'method_id'="SSCOM", MaxCt=MaxCt))
  ##
  res.list<-list(list(Results=Results))
  ## not elegant naming of results list below
  ## but need to have name somewhere as opposed
  ## to just the index to make sure we have the correct stock
  names(res.list)<-StockName
  return(res.list)
}
