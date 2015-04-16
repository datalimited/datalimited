######################################
#RunRegressions--------------------------------------------------
#This code runs each of the desired regressions on supplied regression formatted data
######################################

RunRegressions<- function(Data,RegList,FigureName)
{
  
  #   Data<- RamData
  #   
  #   RegList<- Regressions
  #   
  #   FigureName<- 'eh'
  
  RegNames<- names(RegList)
  
  Models<- list()
  
  for (m in 1:length(RegNames))
  {
    
    ModelVars<- eval(parse(text=paste('RegList$',RegNames[m],sep='')))
    
    #   WhereVars<- colnames(Data) %in% c('IdOrig',ModelVars)
    #   
    WhereVars<- colnames(Data) %in% c(ModelVars)
    #   
    fmla <- as.formula(paste("LogBvBmsy ~ ", paste( ModelVars[ModelVars!='LogBvBmsy' & ModelVars!='IdOrig'], collapse= "+")))
        
    RegData<- Data[Data$ExtendedTime==F,]
    
    UsedNames<- paste(RegData$IdOrig,RegData$Year,sep='-')
    
    AllNames<- paste(Data$IdOrig,Data$Year,sep='-')
    
    TempReg<- lm(fmla,data=RegData)
    
    #   TempReg<- plm(fmla , data=Data[,WhereVars],model="pooling",index=c('IdOrig','Year'),na.action='na.omit')
       
    Model<- RegNames[m]  
        
    Omitted<- names(na.action(TempReg))
    
    HasPrediction<- UsedNames[(rownames(RegData) %in% Omitted)==F]
    
    WherePrediction<- AllNames %in% HasPrediction

    eval(parse(text=paste('Data$',Model,'Marker[WherePrediction]<- TRUE' ,sep='')))
    
    eval(parse(text=paste('Data$',Model,'Prediction[WherePrediction]<- TempReg$fitted.values' ,sep='')))

    # z=predict(TempReg,Data[,WhereVars])
    # 
    # r<- TempReg$residuals
    # 
    # s<- TempReg$model$LogBvBmsy
    
    pdf(file=paste(FigureFolder,FigureName,' ',RegNames[m],' Observed vs Predicted.pdf',sep=''))
    plot(TempReg$fitted.values+TempReg$residuals,TempReg$fitted.values,xlab='Observed Log B/Bmsy',
         ylab='Predicted Log B/Bmsy',main=paste('Model ',RegNames[m],sep=''))
    abline(0,1)
    dev.off()
    
    #   pdf(file=paste(FigureFolder,RegNames[m],' Leverage Plots.pdf',sep=''))
    #   leveragePlots(TempReg,ask=F)
    #   dev.off()
    #   
    show(summary(TempReg))
    
    eval(parse(text=paste('Models$',RegNames[m],'=TempReg',sep='')))
    
  }
  
  return(list(Models=Models,RamData=Data))
} #Close function