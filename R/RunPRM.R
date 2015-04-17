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
load('extra/Ram Data.Rdata')
load('extra/FAO Data.Rdata')

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

sapply(list.files(pattern="[.]R$", path="extra/PRM Functions", full.names=TRUE), source) #load required functions

RamData<- RamData[RamData$Year<=BaselineYear,]

RegressionResults<- RunRegressions(RamData,Regressions,'Real Stocks') # run regressions using all RAM stock values

RealModels<- RegressionResults$Models # use regression models estimated from assessment BvBmsy values only

RealModelFactorLevels<- NULL

Models<- names(Regressions)

TempOmitted<- NULL
for (m in 1:length(names(Regressions)))
{
  Model<- names(Regressions)[m]
  eval(parse(text=paste('RealModelFactorLevels$',Model,'<- RealModels$',Model,'$xlevels$SpeciesCatName',sep='')))
}


RealModelSdevs<- CreateSdevBins(RealModels,RamData,TransbiasBin)


RealModelSdevs<- CreateSdevBins(RealModels,RegressionResults$RamData,TransbiasBin)

WhereFaoNeis<- (grepl('nei',FaoData$CommName) | grepl('spp',FaoData$SciName)) & grepl('not identified',FaoData$SpeciesCatName)==F #Find unassessed NEIs

WhereFaoMarineFish<- grepl('not identified',FaoData$SpeciesCatName)

FaoSpeciesLevel<- FaoData[WhereFaoNeis==F & WhereFaoMarineFish==F ,] #Fao stocks named to the species level

TempLevel<- NULL

TempModel<- NULL

AllPossible<- unique(data.frame(I(FaoData$SpeciesCatName),I(FaoData$SpeciesCat)))

colnames(AllPossible)<- c('SpeciesCatNames','SpeciesCat')

RamPossibleCats<- unique(RamData$SpeciesCatName)

FaoSpeciesPossibleCats<- unique(FaoSpeciesLevel$SpeciesCatName)

Models<- Models[Models!='M7']

for (m in 1:length(Models)) #Apply models to species level fisheries
{

  TempModelName<- Models[m]

  eval(parse(text=paste('TempLevel<- RealModelFactorLevels$',TempModelName,sep='')))

  eval(parse(text=paste('TempModel<- RealModels$',TempModelName,sep='')))

  ProxyCats<- AssignNearestSpeciesCategory(FaoSpeciesLevel,TempLevel,AllPossible)

  Predictions<- predict(TempModel,ProxyCats$Data)

  eval(parse(text=paste('FaoSpeciesLevel$',TempModelName,'Prediction<- Predictions',sep='')))
}


BiomassData<- FaoSpeciesLevel #Only store fisheries that have some form of biomass estimates

BiomassData$BvBmsy<- NA

BiomassData$BestModel<- NA

BiomassData$LogBvBmsy<- NA

BiomassData$FvFmsy<- NA

BiomassData$MSY<- NA

BiomassColumns<- (grepl('BvBmsy$',colnames(BiomassData)) | grepl('Prediction',colnames(BiomassData))) & grepl('LogBvBmsy',colnames(BiomassData))==F

AvailableBio<- (BiomassData[,BiomassColumns])

HasSomething<-  rowSums(is.na(AvailableBio))<dim(AvailableBio)[2]

BiomassData<- BiomassData[HasSomething,]

AvailableBio<- AvailableBio[HasSomething,]

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

Status<- AnalyzeFisheries(BiomassData,'JackStat','Year',min(BiomassData$Year):max(BiomassData$Year),RealModelSdevs,NeiModelSdevs,TransbiasBin,TransbiasIterations)

BiomassData[,c('PrmB')]<- Status$Data$BvBmsy

