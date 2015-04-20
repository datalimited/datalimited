### Run Panel Regression Model (PRM) in the manner of Costello et al. 2012 (and Costello et al. In Prep)
## Written by Dan Ovando, April 16 2015
##########
# This code is actually a Jackknifing routine developed for Costello et al. In Prep.
#The script takes properly formatted RAM data "Ram Data", and
# 1. Sequentially for each fishery...
# 2. Estiamtes six PRM models in the manner of Costello et al. 2012 omitting a given gishery
# 3. Uses the PRM to estimate B/Bmsy for the omitted fishery
# 4. Finds the "best" (defined as most data rich) model that can be used for that fishery, and stores that result
#
###########
rm(list=ls())
load('Ram Data.Rdata')
NumCPUs<- 1 #Can be run in parallel if desired, not really useful for jackknifing one at a time
BaselineYear<- 2012 #Leave this as 2012
FigureFolder<- paste('extra/PRM Figures/')
dir.create(FigureFolder,recursive=T)

DependentVariable<- 'BvBmsy' #Dependent variable in regression

IsLog<- TRUE #Should dependent variable be logged?

DependentName<- DependentVariable

DependentName<- if (IsLog==T){paste('Log',DependentVariable,sep='')}

CatchLags<- 4 #Number of years of lagged catch to create for regression, leave this alone for now

LifeHistoryVars<- c('MaxLength','AgeMat','VonBertK','Temp') #Life history variables to include for potential regression

IdVar<- 'IdOrig' #Id variable to use in regressions

CatchVariables<- c('YearsBack','ScaledCatch',paste('ScaledCatch',1:CatchLags,'Back',sep=''),'MaxCatch','TimeToMaxCatch','InitialScaledCatchSlope'
                   ,'MeanScaledCatch','CatchToRollingMax')

Regressions<- list(M1=c(DependentName,CatchVariables,LifeHistoryVars,'SpeciesCatName'),M2=c(DependentName,CatchVariables,'MaxLength','AgeMat','VonBertK','SpeciesCatName'),
                   M3=c(DependentName,CatchVariables,'MaxLength','VonBertK','SpeciesCatName'),M4=c(DependentName,CatchVariables,'VonBertK','SpeciesCatName'),M6=c(DependentName,CatchVariables,'SpeciesCatName'))

ModelNames<- names(Regressions) #Create columns to store the results of each PRM


TransbiasBin<- 0.9

library(plyr)
library(lattice)
library(rfishbase)
library(stringr)
library(RCurl)
library(XML)
library(MASS)
library(prettyR)
library(zoo)
library(proftools)
library(snowfall)
library(parallel)
# library(shiny)
library(ggplot2)
library(gridExtra)
library(reshape2)

sapply(list.files(pattern="[.]R$", path="PRM Functions", full.names=TRUE), source) #load required functions

RamData<- RamData[RamData$Year<=BaselineYear,]

Regions<- unique(RamData$Country)

RamIds<- unique(RamData$IdOrig)

JackStore<- as.data.frame(matrix(NA,nrow=0,ncol=14))

colnames(JackStore)<- c('Assessid','Year','Country','Catch','RamB','RamF','RamMSY','PrmB')

TransbiasIterations<- 1000

RamIds<- unique(RamData$IdOrig)

for (r in 1:length(RamIds))
{

  Omit<- RamData[RamData$IdOrig%in%RamIds[r],] #Find the fishery to omit

  if (sum(is.na(Omit$Catch))==0) # If there's any catch to work with
  {

    Omit$CatchToRollingMax[is.na(Omit$CatchToRollingMax)]<- 0

    Jacked<- RamData[!(RamData$IdOrig%in%RamIds[r]),]

    TempJack<- as.data.frame(matrix(NA,nrow=dim(Omit)[1],ncol=14))

    colnames(TempJack)<- c('Assessid','Year','Country','Catch','RamB','RamF','RamMSY','PrmB')

    TempJack[,c('Assessid','Year','Country','Catch','RamB','RamF','RamMSY')]<- Omit[,c('IdOrig','Year','Country','Catch','BvBmsy','FvFmsy','MSY')]

    JackModel<- RunRegressions(Jacked,Regressions,'Real Stocks')

    RealModelFactorLevels<- NULL

    Models<- names(Regressions)

    TempOmitted<- NULL

    for (m in 1:length(names(Regressions)))
    {
      Model<- names(Regressions)[m]
      eval(parse(text=paste('RealModelFactorLevels$',Model,'<- JackModel$Models$',Model,'$xlevels$SpeciesCatName',sep='')))
    }

    Jacked<- InsertFisheryPredictions(Jacked,JackModel) #Add fishery predictions back into main dataframe

    RealModelSdevs<- CreateSdevBins(JackModel$Models,Jacked,TransbiasBin)

    ## This section of code assigns the nearest available species category to any stocks of a
    ## of a species category not encompassed by the regression
    AllPossible<- unique(data.frame(I(Jacked$SpeciesCatName),I(Jacked$SpeciesCat)))

    colnames(AllPossible)<- c('SpeciesCatNames','SpeciesCat')

    RamPossibleCats<- unique(RamData$SpeciesCatName)

    Models<- Models[Models!='M7']

    for (m in 1:length(Models)) #Apply models to species level fisheries
    {

      TempModelName<- Models[m]

      eval(parse(text=paste('TempLevel<- RealModelFactorLevels$',TempModelName,sep='')))

      eval(parse(text=paste('TempModel<- JackModel$Models$',TempModelName,sep='')))

      ProxyCats<- AssignNearestSpeciesCategory(Omit,TempLevel,AllPossible)

      Predictions<- predict(TempModel,ProxyCats$Data)

      eval(parse(text=paste('Omit$',TempModelName,'Prediction<- Predictions',sep='')))
    }

    BiomassData<- Omit #Only store fisheries that have some form of biomass estimates

    BiomassData$BvBmsy<- NA

    BiomassData$BestModel<- NA

    BiomassData$LogBvBmsy<- NA

    BiomassData$FvFmsy<- NA

    BiomassData$MSY<- NA

    BiomassColumns<- (grepl('BvBmsy$',colnames(Omit)) | grepl('Prediction',colnames(Omit))) & grepl('LogBvBmsy',colnames(Omit))==F

    AvailableBio<- (BiomassData[,BiomassColumns])

    HasSomething<-  rowSums(is.na(AvailableBio))<dim(AvailableBio)[2]

    BiomassData<- BiomassData[HasSomething,]

    AvailableBio<- AvailableBio[HasSomething,]

    TempJack<- TempJack[HasSomething,]

    AvailableBioMarker<- matrix(rep((1:dim(AvailableBio)[2]),dim(AvailableBio)[1]), dim(AvailableBio)[1],dim(AvailableBio)[2],byrow=TRUE)

    AvailableBioMarker<- AvailableBioMarker*(is.na(AvailableBio)==F)

    AvailableBioMarker[AvailableBioMarker==0]<- NA

    BestModel<- apply(AvailableBioMarker,1,min,na.rm=T)

    BestBio<- NULL
    for (b in 1:dim(AvailableBio)[1])
    {
      BestBio[b]<- AvailableBio[b,BestModel[b]]
    }

    BestBio[BestModel==1]<- log(BestBio[BestModel==1])

    BestModelnames<- c('RAM',ModelNames)

    BestModelNames<- BestModelnames[sort(unique(BestModel))]

    BestModel<- as.factor((BestModel))

    levels(BestModel)<- BestModelNames

    BiomassData$BestModel<- BestModel

    BiomassData$BvBmsy<- BestBio

    BiomassData$IdLevel<- 'Species'

    BiomassData$Dbase<- 'FAO'

    OmitStatus<- AnalyzeFisheries(BiomassData,'JackStat','Year',min(BiomassData$Year):max(BiomassData$Year),RealModelSdevs,NeiModelSdevs,TransbiasBin,TransbiasIterations)

    TempJack[,c('PrmB')]<- OmitStatus$Data$BvBmsy

    show(paste(100*(r/length(RamIds)),' % Done with JackKnife',sep=''))

    JackStore<- rbind(JackStore,TempJack)
  } #Close if all catch loop
} #Close stock loop


Prm<- cbind(JackStore[,c('Assessid','Year','Country','Catch','RamB','PrmB','RamMSY','CmsyMSY','RamF','CmsyF')],'PRM')

colnames(Prm)<- c('Id','Year','Country','Catch','RamBvBmsy','ModelBvBmsy','RamMSY','CmsyMSY','RamF','CmsyF','Model')

PlotJack<- Prm

SpeciesInfo<- RamData[,c('IdOrig','SpeciesCatName','MaxLength','AgeMat','VonBertK')]

colnames(SpeciesInfo)<- c('Id','SpeciesCatName','MaxLength','AgeMat','VonBertK')

PlotJack<- join(PlotJack, SpeciesInfo,by='Id',match='first')

PlotJack$ProportionalBError<- 100*((PlotJack$ModelBvBmsy-PlotJack$RamBvBmsy)/PlotJack$RamBvBmsy)


pdf(file=paste(FigureFolder,'Proportional error in BvBmsy by time.pdf',sep=''))

B_Error_By_Time<- ( ggplot(data=subset(PlotJack,Year>1985 & Model=='Cmsy'),aes(factor(Year),ProportionalBError))+
                      geom_boxplot(outlier.shape=NA,fill='steelblue2')+
                      coord_cartesian(ylim=c(-200,200))+
                      #   coord_cartesian(ylim = range(boxplot(PlotJack$ProportionalBError, plot=FALSE)$stats)*c(.8, 1.5))+
                      geom_abline(intercept=0,slope=0)+
                      ylab('Proportional Error (%) in B/Bmsy')+
                      xlab('Time')+
                      scale_x_discrete(breaks=as.character(seq(from=1986,to=2012,by=8))))

print(B_Error_By_Time)

dev.off()


