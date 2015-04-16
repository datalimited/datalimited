######################################
#InsertFisheryPredictions--------------------------------------------------
# This code adds fisheries that were able to obtain predictions of B/Bmsy using PRM to primary data frame 
######################################

InsertFisheryPredictions<- function(Data,Models)
{

#   Data<- RamData
#   
#    Models<- RealModels
# # 
   ModelNames<- names(Models)
  
  for (m in 1:length(ModelNames))
  {
    
  Model<- ModelNames[m]  
    
  eval(parse(text=paste('TempOmitted<- na.action(Models$',Model,')',sep='')))
  
  KeptEntries<- ((1:dim(Data)[1]) %in% TempOmitted)==F #Identify fisheries used in each model
  
#   eval(parse(text=paste('Data$',Model,'Marker<- FALSE' ,sep='')))
  
  eval(parse(text=paste('Data$',Model,'Marker[KeptEntries]<- TRUE' ,sep='')))
  
  eval(parse(text=paste('Data$',Model,'Prediction[KeptEntries]<- Models$',Model,'$fitted.values' ,sep='')))
  } #Close model loop
  
  return(Data)
} #Close function 