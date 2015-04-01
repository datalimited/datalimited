#################COM SIR FUNCTIONS #################################
#  Catch only model implemented by C.V.Minte-Vera                                             #
#  Version June  11, 2013                                                                                 #
#  Vasconcellos and Cochrane 2005                                                                    #
#  Biomass dynamic (Schaeffer) and harvest dynamic model (logistic)                    #
#  assumed the initial harvest rate equals Catch first year / Biomass first year          #
#  Biomass first year = carrying capacity                                                              #
#  lognormal likehood for catch data w/ observation error CV = 0.4                          #
#  Estimation using Bayesian Sampling Importance Ressampling                             #
#  Joint prior is the important function, i.e. resampling proportional to the likelihood   #
# changes: cap on harvest rate
# prior on x<-runif(Nsim,0.000001,1)# prior for x
#################   <><   <><   <><   ####################################

CheckData<-function(MyData){
for (i in 1:dim(MyData)[2]) cat(colnames(MyData[i])," ",class(MyData[,i])," ", table(MyData[,i]),"\n")
}

#This generates prior values  INCLUDE RESILIENCE AS
# prior for r and  logKl~U(log(max(catch)) , log(100*max(catch))) as prior for K
# following Catch-MSY assumptions

#change prior for x
Priors<-function(Nsim,logK=T,MyPar=TruePar,CV=0.4,NormK=F,Normr=F,
			Norma=F,Normx=F,Catch,minK=max(Catch), maxK=100*max(Catch),start.r=c(0.2,1),LogisticModel=T,Obs=F)
{
	if(logK)
	#&& NormK==F) #uniform prior on the log(K)
    {
	        K<- runif(Nsim,log(minK),log(maxK))
	        K<- exp(K)
   	}
	else K<-runif(Nsim,minK,maxK)# prior on arithmetic scale for K

   	if(logK==F && NormK==T) #normal prior on the arithmetic scale
       	{K<-rnorm(Nsim,MyPar$K,CV*MyPar$K)}
	if(Normr) r<-rnorm(Nsim,MyPar$r,CV*MyPar$r)
	else r<-runif(Nsim,start.r[1],start.r[2])# prior for r
	if(Normx) x<-rnorm(Nsim,MyPar$x,CV*MyPar$x)
	#else x<-runif(Nsim,0,10)# prior for x
	#else x<-runif(Nsim,0,0.5)# prior for x
	else x<-runif(Nsim,0.000001,1)# prior for x

	if(Norma) a<-rnorm(Nsim,MyPar$a,CV*MyPar$a)
	else a<-runif(Nsim,0,1)# prior for a
	#sigma<-rnorm(Nsim,0.001,0.0001)# normal informative prior for sigma
	#h<-runif(Nsim,0,0)# normal informative prior for initial harvest rate
    h<-rep(0,Nsim)# h=0, start from unexploited state
	z<-rep(1,Nsim)# z=1 is the Schaeffer model

	N1<-(1-h)*K #initial population size
    Like<-rep(1,Nsim) #Like=1 if the model does not crashes, Like=0 if the model crashes

	#Set up the numbers and project to the start of the first "real" year
	#(i.e. N1 to Nm), all simulations at once
	#initial conditions
	predbio<- N1
	predcatch<-Catch[1] #the first year of catches in assumed known
	predprop<-predcatch/predbio
	inipredprop<-predprop
	m<-length(Catch)
	 for (t in 1:m-1)
     {
	    #effort dynamics
	   	if (LogisticModel) #logistic model with Bt-1
		     predprop<- predprop*(1+x*((predbio/(a*K)-1)))

	   		#predprop = predprop*(1+x*((predbio/(K*a)) -1))
	   	else #linear model
	   		predprop = predprop + x*inipredprop;
	   	#biomass dynamics
	   	if(Obs) predbio = predbio + ( r * predbio * ( 1 -( predbio / K ))) - Catch[t]
	   	else predbio = predbio + ( r * predbio * ( 1 -( predbio / K ))) - predcatch
	   	predcatch = predbio*predprop
        for(i in 1:Nsim) #this is an element operation ## change 07 Jan 2013 "|| is.na(Like[i])"
			if ((Like[i]!=0 & (predbio[i]<= 0 || predprop[i]< 0 )) || is.na(Like[i]) ) Like[i]<- 0


	  }


	#in some of the years the population crashes, we know this is impossible, because the population
	# is still there (because of the catches
	#Output
	Params<-NULL
	Params$N1<-N1
	Params$K<-K
	Params$r<-r
	Params$z<-z
	Params$a<-a
	Params$x<-x
	Params$h<-h
	Params$B<- predbio
	Params$Prop<-predprop
    Params$Like<-Like
	return(as.data.frame(Params))
}




#x3<-Priors(MyPar=SWO.South.Start,Catch=t(ct),LogisticModel=T,Nsim=1000,logK=T,NormK=F,Normr=F,Obs=F,start.r=c(0.01,1.8))


###########18 March 2013##############

MYsumpars<- function(ParVals,Npost=1000,Catch,Plot=F)
{
	#Quick version of SIR #ressampling
	#print(ParVals)
	Nsim<- length(ParVals$N1)
	#print(paste("simulations in ",Nsim))
	#print(paste("simulations out ",Npost))

	if(Plot)
	{
		###new
		 tiff(file="hist.tif",width=4,height=4,units="in",res=200)

			#par(mfrow=c(2,4),mgp=c(1.5,0.5,0), mar=c(0.5,0.5,0.5,0.5), cex=0.4,tck=-0.02)
		par(mfrow=c(2,4), mar=c(3.5,1,1,1), cex=0.4,tck=-0.02)

				#plot priors
			hist(ParVals$K,xlab="K",main="K no crash only")
			hist(ParVals$r,xlab="r",main=" r ")
			hist(ParVals$a,xlab="a",main=" a ")
			hist(ParVals$x,xlab="x",main=" x  ")
	}
	# if NA the prob will be zero 20 Jan 2013
	#print("how many probabilities are NAs?")
	#print(table(is.na(ParVals$Like))) #debug
	NA.Prob<-sum(as.numeric(is.na(ParVals$Like)))
	ind1<-which(is.na(ParVals$Like))
	ParVals$Like[ind1]<-rep(0,length(ind1))

	Ind<- sample(Nsim,Npost,replace=T,prob=ParVals$Like) #sample with replacement (just the index)
	    # this is an R function will sample proportional to the likelihood
		#'sample' takes a sample of the specified size from the elements of
	    #'x' using either with or without replacement.
	Ind<-sort(Ind)
	m<-length(ParVals)
	ParVals<-as.data.frame(ParVals)
	#[1] "N1"   "K"    "r"    "z"    "a"    "x"    "h"    "B"    "Prop" "Like"
	Post<- ParVals[Ind,]
	#colnames(Post)<-names(ParVals)
	#plot posteriors
	if(Plot)
	{	hist(Post[,2],xlab="K",main="posteriors")
		hist(Post[,3],xlab="r",main=" ")
		hist(Post[,5],xlab="a",main="  ")
		hist(Post[,6],xlab="x",main="  ")

	dev.off()
	#print(Ind)
	}
	#time series of BoverBmsy
	 #BoverBmsy=apply(MARGIN=1,FUN=effortdyn_orig,X=Post[,1:7],Catch=Catch,LogisticModel=T)
	#
	BoverBmsy=apply(MARGIN=1,FUN=effortdyn_orig,X=Post[,1:7],Catch=Catch,LogisticModel=F,What=1)
	#=xx$BoverBmsy
	#BoverBmsy","biomass","predprop","predcatch","obscatch"
	biomass=apply(MARGIN=1,FUN=effortdyn_orig,X=Post[,1:7],Catch=Catch,LogisticModel=T,What=2)
	predprop=apply(MARGIN=1,FUN=effortdyn_orig,X=Post[,1:7],Catch=Catch,LogisticModel=T,What=3)
	predcatch=apply(MARGIN=1,FUN=effortdyn_orig,X=Post[,1:7],Catch=Catch,LogisticModel=T,What=4)
	residual= Catch-predcatch

	#BoverBmsy<-matrix(NA,ncol=length(Catch),nrow=Npost)
	#for(j in 1:Npost) BoverBmsy[j,]=effortdyn_orig(Post[j,1:7],Catch=Catch,LogisticModel=T)

	#diagnostics
	repeatedSamples<-table(Ind)
	repeatedSamples<- table(repeatedSamples)#number of Ind with 1, 2 or more samples
	#MSD - maximum single density,should be less than 1%
	MSD<- (max(as.numeric(names(repeatedSamples)))/Npost)*100

	##novo
	MIR<- max(Post$Like)/sum(Post$Like) ###maximum importance ratio or maximm importantance weight
	CV.IR<-  ((1/length(Post$Like))*sum((Post$Like)^2)) - ((1/length(Post$Like))*sum(Post$Like))^2
	CV.IR<- sqrt(CV.IR)/((length(Post$Like)^(-0.5))*sum(Post$Like))
	#print("MSD - maximum single density,should be less than 1%")
	#print(repeatedSamples)
	#print("MIR - maximum importance ratio,should be less than 0.04")
	#print(MIR)
	#print("CV.IR - CV importance ratio, should be less than 0.04")
	#print(CV.IR)
	##new lines 06 Jan 2013  create a file with the results
	#write.table(repeatedSamples, "MSD.txt")
	#write.table(MIR, "MIR.txt")
	#write.table(CV.IR, "CVIR.txt")
	###end
	#Raftery and Bao 2010 diagnostics - CHECK
	#(1) Maximum importance weight = MIR
	# (2) variance of the importance weights
	Weights<-Post$Like/sum(Post$Like)
	Var.RW= (Npost*Weights -1)^2 #vector
	Var.RW= 	(1/Npost)*sum(Var.RW) #scalar
	# (3) entropy of the importance weights relative to uniformity
	Entropy= -1* sum(Weights*(log(Weights)/log(Npost)))
	# (4)Expected number of unique points after resampling
	Exp.N= sum((1-(1-Weights)^Npost))
	# Effective sample size ESS
	ESS= 1/(sum(Weights^2))


	diagnostics=c(Nsim,NA.Prob,MIR,CV.IR,MSD,Var.RW,Entropy,Exp.N,ESS)
	names(diagnostics)=c("N.nocrash.priors","NA.Prob","MIR","CV.IR","MSD","Var.RW","Entropy","Exp.N","ESS")
	#return(as.data.frame(Post))
    #return(list(posterior=Post, BoverBmsy=BoverBmsy,Diagno=diagnostics))
    return(list(posterior=Post, BoverBmsy=BoverBmsy,Biomass=biomass,E=predprop,C.hat=predcatch,Res=residual,Diagno=diagnostics))

}

#########################ORIGINAL FUNCTIONS
Estimate<-function(Post,Catch,CV=0.4,LogisticModel=F,NormalL=T)
{
	#Post=yy;CV=0.4; Catch=MyCatch;LogisticModel=F;NormalL=F
	#NormallL=T #for normal likelihood
	# NormalL=F for lognormal likelihood
		#The log normal distribution has density
	#  f(x) = 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))
	#     where mu and sigma are the mean and standard deviation of the
	#     logarithm. The mean is E(X) = exp(mu + 1/2 sigma^2), and the
	#     variance Var(X) = exp(2*mu + sigma^2)*(exp(sigma^2) - 1) and hence
	#     the coefficient of variation is sqrt(exp(sigma^2) - 1) which is
	#     approximately sigma when that is small (e.g., sigma < 1/2).
			#x<-seq(0.001,100,0.001)
			#y<-dlnorm(log(x),meanlog=log(1),sdlog=1,log=FALSE)
			# dlnorm(1) == dnorm(0)
			#z<-dnorm(log(x),mean=0,sdlog=1,log=FALSE)
			#plot(log(x),y,type="l")
			#plot(log(x),z,type="l")

    Nsim<- length(Post$K)
	#print(paste("Number of simulations ",Nsim,sep=""))
	# the parameters come from the resampling part
	N1<- Post$N1
	K<- Post$K
	r<- Post$r
	z<- Post$z
	a<- Post$a
	x<- Post$x
	h<- Post$h
	Like<- Post$Like
	#sigma<-rnorm(Nsim,0.001,0.0001)# normal informative prior for sigma
	#Set up the numbers and project to the start of the first "real" year
	#(i.e. N1 to Nm), all simulations at once
	#initial conditions
	m<-length(Catch)
	Params<-NULL
	logLike<- rep(0,Nsim)
	pen1<-rep(0,Nsim)
	pen2<-rep(0,Nsim)
	temp1<-matrix(0,ncol=2,nrow=Nsim)
	temp2<-matrix(0,ncol=2,nrow=Nsim)
	CumLogLike<-rep(0,Nsim)
	for (t in 1:m)
    {
	    if(t==1)
	    {
		    predbio<- N1
			predcatch<-Catch[1] #first catch is assumed known -- BIG ASSUMPTION ---
			predprop<-predcatch/predbio
			inipredprop<-predprop
		}
	    else
	    {
		    #effort dynamics
		   	if (LogisticModel) #logistic model with Bt-1
		   		temp1<- sapply(predprop*(1+x*((predbio/(a*K)-1))),FUN=posfun_orig,eps=0.00001,simplify=TRUE)
		   		#predprop <- predprop*(1+x*((predbio/(K*a)) -1))
		    else #linear model
		   		temp1<-  sapply((predprop + x*inipredprop),FUN=posfun_orig,eps=0.00001,simplify=TRUE)
		   	predprop<-pmin(temp1[1,],0.99) #maximum harvest rate is alwasy 0.99
		   	pen1<- pen1 + temp1[2,]

		   	#biomass dynamics
		    temp2 <- sapply(predbio + ( r * predbio * ( 1 -( predbio / K )))- predcatch,FUN=posfun_orig,eps=0.00001,simplify=TRUE)
		    predbio<-temp2[1,]
            pen2<- pen2 + temp2[2,]
		   	predcatch <-  predbio*predprop
		   	#print(c(temp1,temp2))
   		}
	   	# assumption of CV=0.4 in Vasconcellos and Cochrane
	   	if(NormalL)
	   	{
		   	Mysd<-CV*predcatch
	   		#logLike<- logLike + -0.5*((Catch[t]-predcatch)/Mysd)^2 ##cummulative normal log likelihood
	   		#print(c(t,Catch[t],predcatch,logLike))
			logLike<-dnorm(Catch[t],mean=predcatch,sd=Mysd,log=TRUE)##normal Log likelihood
	   		CumLogLike<- CumLogLike + logLike ## cummulative normal (complete) loglikelihood
	   		#print(c(logLike,CumLogLike))
	   	}
   		else
   		{
	   	 	#logLike<- logLike - log(Catch[t]) -0.5*(log(Catch[t+1])-log(predcatch)/CV)^2
	   		#print(c(t,Catch[t],predcatch,logLike))
			logLike<- dlnorm(Catch[t],meanlog=log(predcatch), sdlog=CV, log=TRUE)
			CumLogLike<- CumLogLike + logLike
	   		#print(c(logLike,CumLogLike))


		}
		# print(c(predbio,predprop,predcatch,Catch[t]))
	  # print(c(predbio))
	}
	#print(c(CumLogLike,pen1,pen2))
    logLike<- rep(0,Nsim) # set vector to 0
    logLike<- CumLogLike - pen1 - pen2 #the penalties turn out too big
    #logLike<- CumLogLike
    #logLike<- logLike - max(logLike) #re-center with logLike for MLE = 0
    TotLike<- sum(logLike,na.rm=TRUE)
    Like<- exp(logLike)/sum(exp(logLike)) #re-normalize so it will add up to 1
	#Like<-(exp(logLike)/exp(TotLike))*100 #re-normalize the likelihood used for re-sampling in SIR
 	#in some of the years the population crashes, we know this is impossible, because the population
	#is still there (because of the catches
	#Output
	Params$N1<-N1
	Params$K<-K
	Params$r<-r
	Params$z<-z
	Params$a<-a
	Params$x<-x
	Params$h<-h
  	Params$B<-predbio #biomass in the latest year
  	Params$Prop<-predprop # harvest rate in the latest year
	Params$LogLike<- CumLogLike #Loglikelood used for DIC
 	Params$Like<-Like #re - normalized likelihood for SIR
    return(Params)
}


posfun_orig<- function(x,eps)
{	# returns a newvalue to replace the old and a penalty for the likelihood
	#this adds a penalty


    if (x>=eps)
        return(c(x,0))
    else
    {
	   pen<-.01*(x-eps)^2
       y<- 1.0-x/eps
       newvalue<-eps*(1.0/(1.0+y+y*y))
       return(c(newvalue,pen))
    }
}



effortdyn_orig<-function(TrueP,Catch,LogisticModel,What)
{
	#What: what to return
	# What=1   "BoverBmsy",
	# What=2 "biomass",
	# What=3  "predprop",
	# What=4 "predcatch",
	# What=5 "obscatch"

	#True Parameters
	hrate=TrueP['h']
	 K=TrueP['K']
	 r=TrueP['r']
	 x=TrueP['x']
	 a=TrueP['a']
	 m=length(Catch)
	 predbio=predprop=predcatch= rep(0,m)

	 #print(TrueP)
	 #initial conditions
	 predbio[1]=(1.0-hrate)*K
	 predprop[1]=Catch[1]/predbio[1]
	 predcatch[1]= Catch[1]

	 for (t in 2:m)
     {
		#biomass dynamics
	   	predbio[t] = predbio[t-1]+ (r*predbio[t-1]*(1-(predbio[t-1]/K))) - predcatch[t-1]
	   	#effort dynamics
	   	if (LogisticModel) #logistic model with Bt-1
	   		predprop[t] = predprop[t-1]*(1+x*((predbio[t-1]/(K*a)) -1))
	   	else #linear model
	   		predprop[t] = predprop[t-1] + x*predprop[1];
	   	predcatch[t] = predbio[t]*predprop[t]
     }

     BMSY= K/2.0
     BoverBmsy=predbio/BMSY

     xx<-cbind(BoverBmsy,predbio,predprop,predcatch,Catch)
     names(xx)<-c("BoverBmsy","biomass","predprop","predcatch","obscatch")
     #return(BoverBmsy)
     return((xx[,What]))

 }




############ NEW

DoProject<-function(Simul=T,MyFile="Figura1",Myseed=128,TruePar,MyData,MyPri,EstLogisticM=F,LogisticModel=F,MyCV=0.4,NormalL=F,Nsim=200,Npost=100,MyYLim = c(-1, 10),...)
{
	#This functions calls others that implement SIR to estimate the posteriors distributions
	set.seed(Myseed)
	names(TruePar)<-c( "N1",   "K",    "r",    "z",    "a",    "x",    "h")
	TruePar<-as.data.frame(t(TruePar))
	MyCatch<- MyData$catch
	ww<-Priors(Nsim=Nsim,MyPar=TruePar,Catch=MyCatch,...)#sample from the priors & compute Likelihood L=0 if model crashes, L=1 o/w
	#print(paste("simulations in ",length(ww[,1])))
	yy<-ww[ww$Like>0,]
	#print(paste("simulations out ",length(yy$Like)))
	hh<-Estimate(Post=yy,CV=MyCV,Catch=MyCatch,LogisticModel=EstLogisticM,NormalL=NormalL)#compute Likelihood for those models that do not crash
	zz<-MYsumpars(Npost=Npost,ParVals=hh,Catch=MyCatch)# resample proportional to likelihood
	#write.table(zz,"post.txt")
	#write.table(yy,"posmodelprior.txt")
	#write.table(ww,"priors.txt")

  zz
}


