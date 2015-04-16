######################################
#Retransformation Bias--------------------------------------------------
# This code applies the retransformation bias to any fisheries that need it
######################################



TransBias<- function(Data,SdevBins,BinBreak,J)
{
  
#    Data<-  subset(BiomassData,(IdLevel=='Neis' | IdLevel=='Unidentified') & Dbase=='FAO')
#    
#   SdevBins<- NeiModelSdevs
#   
#   BinBreak<- TransbiasBin
#   
#   J<- 100
  
  core<- Data[,c('IdOrig','BvBmsy','BestModel','Year')]
  
  colnames(core)<- c('id','bvb','mod','years')
  
  ModelNames<- names(SdevBins)
  
  #   stdbin$high[is.na(stdbin$high)]<- 0
  #   stdbin$low[is.na(stdbin$low)]<- 0
  #   originalstd<- stdbin
  
  um<- sort(unique(core$mod)) #unique models in core
  
  core$bvb<- exp(core$bvb) #exponentiate raw bvbmsy values 
  
  years<- sort(unique(core$years)) #unique years
  
  fs<- sort(unique(core$id)) #unique fishery ids
  
  bt<- matrix(NA,nrow=dim(core)[1],ncol=1) #blank matrix for results
  
  bt[core$bvb>=BinBreak,1]<- 1 #mark bin locations
  
  bt[core$bvb<BinBreak,1]<- 0 #mark bin locations
  
  core$bin<- bt
    
  indtemp<- matrix(NA,nrow=dim(core)[1],ncol=9) #Matrix for individual fishery results
  
  indtemp<- as.data.frame(indtemp)
  
  colnames(indtemp)=c('id','years','raw','mean','top','bot','iq75','iq25','logsd')
  
  c<- matrix(NA,nrow=length(years),ncol=6)
  
  c<- as.data.frame(c)	
  
  colnames(c)<- c('year','median','top','bottom','iq75','iq25')
  
  bin<- matrix(NA,nrow=3,ncol=3)
  
  bin<- as.data.frame(bin)
  
  colnames(bin)=c('O','F','U')
  
  rownames(bin)=c('med','bot','top')
  
  # bin<- matrix(NA,nrow=3,ncol=4)
  # bin<- as.data.frame(bin)
  # colnames(bin)=c('<4','4to8','8to12','>12')
  # rownames(bin)=c('med','bot','top')
  
  collapsed<- matrix(NA,nrow=length(years),ncol=4)
  
  colnames(collapsed)=c('year','med','bot','top')
  
  collapsed<- as.data.frame(collapsed)

  medians<- matrix(NA,nrow=length(years),ncol=J)
  
#   colnames(medians)=c('year','med','bot','top')
  
  collapsed<- as.data.frame(collapsed)
  
  under<- matrix(NA,nrow=length(years),ncol=4)
 
  colnames(under)=c('year','med','bot','top')
  
  under<- as.data.frame(under)
  
  over<- matrix(NA,nrow=length(years),ncol=4)
  
  colnames(over)=c('year','med','bot','top')
  
  over<- as.data.frame(over)
  
  PossibleModels<- names(SdevBins)
    
  
  # Correct retransformation bias -------------------------------------------
    
  DistStorage<- as.data.frame(matrix(NA,nrow=0,ncol=J+2))
  
  for (y in 1:length(years)) #loop over years
  {
    
    where<- core$years==years[y] #find year locations
    
    tcore<- core[where,] #temporary core
    
    um<- unique(tcore$mod)
    
    jstore<- matrix(NA,nrow=dim(tcore)[1],ncol=J) #blank matrix for bootstrap
    
    
    
    for (m in 1:length(um)) #loop over models
    {
      
      whereh<- tcore$mod==um[m] & tcore$bin==1 & is.na(tcore$bvb)==F	#find fisheries in high bin		
      wherel<- tcore$mod==um[m] & tcore$bin==0 & is.na(tcore$bvb)==F #find fisheries in low bin	
      
      WhichModel<- c(1:length(PossibleModels))[names(SdevBins)==um[m]]
      
      stdbin<- SdevBins[[WhichModel]]
      
      stdbin$HighBvBmsy[is.na( stdbin$HighBvBmsy)]<- mean(stdbin$HighBvBmsy,na.rm=T)
      
      stdbin$LowBvBmsy[is.na( stdbin$LowBvBmsy)]<- mean(stdbin$LowBvBmsy,na.rm=T)
      
      ### Series of steps to deal with some time issues ###
      # if (years[y]<=max(stdbin$years) & years[y]>=1950)
      if (years[y]<=max(stdbin$Year) & years[y]>=min(stdbin$Year))
      {
        wheres<- stdbin$Year==years[y] 
        
        if (sum(wheres)>0)
        {
        
        highstd<- stdbin$HighBvBmsy[wheres]
        
        lowstd<- stdbin$LowBvBmsy[wheres]		
        }
        else
        {
          highstd<- mean(stdbin$HighBvBmsy,na.rm=T)
          
          lowstd<- mean(stdbin$LowBvBmsy,na.rm=T)
          
        }
      }
      if (years[y]<min(stdbin$Year))
      {
        wheres<- stdbin$Year==min(stdbin$Year) 
        if (sum(wheres)>0)
        {
          
        highstd<- stdbin$HighBvBmsy[wheres]
        
        lowstd<- stdbin$LowBvBmsy[wheres]  	
        }
        else
        {
          highstd<- mean(stdbin$HighBvBmsy,na.rm=T)
          
          lowstd<- mean(stdbin$LowBvBmsy,na.rm=T)
          
        }
      }
      if (years[y]>max(stdbin$Year))
      {
        wheres<- stdbin$Year==max(stdbin$Year)
            
        if (sum(wheres)>0)
        {
          
        highstd<- stdbin$HighBvBmsy[wheres]
        
        lowstd<- stdbin$LowBvBmsy[wheres] 	
        }
        else
        {
          highstd<- mean(stdbin$HighBvBmsy,na.rm=T)
          
          lowstd<- mean(stdbin$LowBvBmsy,na.rm=T)
          
        }
        
      }
      

      ### Apply error terms ###
      
      jstd<- rnorm(sum(whereh,na.rm=T)*J,mean=0,sd=highstd)
      dim(jstd)<- c(sum(whereh,na.rm=T),J)
      #       jstore[whereh,]<- rep(log(tcore$bvb[whereh]),J) +highstd*jstd
      jstore[whereh,]<- rep(log(tcore$bvb[whereh]),J) + jstd
      
      jstd<- rnorm(sum(wherel,na.rm=T)*J,mean=0,sd=lowstd)
      dim(jstd)<- c(sum(wherel,na.rm=T),J)
      #       jstore[wherel,]<- rep(log(tcore$bvb[wherel]),J) +lowstd*jstd
      jstore[wherel,]<- rep(log(tcore$bvb[wherel]),J) + jstd
    } #close loop over models
    
    
    DistStorage<- rbind(DistStorage, data.frame(years[y],as.character(core$id[where]),jstore))
    
    ### Store Results ###
    
    ### % collapsed###
    ctemp<- exp(jstore)<=.2
    ctemp<- colSums(ctemp,na.rm=T)/sum(where)
    csort<- sort(ctemp)
    ctop<- csort[ceiling(0.975*length(csort))]
    cbot<- csort[ceiling(0.025*length(csort))]
    cmed<- median(csort,na.rm=T)
    
    ### % B > Bmsy ###
    otemp<- exp(jstore)<1
    otemp<- colSums(otemp,na.rm=T)/sum(where)
    osort<- sort(otemp)
    otop<- osort[ceiling(0.975*length(osort))]
    obot<- osort[ceiling(0.025*length(osort))]
    omed<- median(osort,na.rm=T)
    
    ### % B<Bmsy ###
    utemp<- exp(jstore)>=1
    utemp<- colSums(utemp,na.rm=T)/sum(where)
    usort<- sort(utemp)
    utop<- usort[ceiling(0.975*length(usort))]
    ubot<- usort[ceiling(0.025*length(usort))]
    umed<- median(usort,na.rm=T)
    
    
    ### % in different bins ###
    b1temp<- exp(jstore)<=0.8
    # b1temp<- exp(jstore)<=0.5
    
    b1temp<- colSums(b1temp,na.rm=T)/sum(where)
    b1sort<- sort(b1temp)
    b1top<- b1sort[ceiling(0.975*length(b1sort))]
    b1bot<- b1sort[ceiling(0.025*length(b1sort))]
    b1med<- median(b1sort,na.rm=T)
    
    # b2temp<- exp(jstore)<=1.5 & exp(jstore)>0.5
    
    b2temp<- exp(jstore)<=1.2 & exp(jstore)>0.8
    b2temp<- colSums(b2temp,na.rm=T)/sum(where)
    b2sort<- sort(b2temp)
    b2top<- b2sort[ceiling(0.975*length(b2sort))]
    b2bot<- b2sort[ceiling(0.025*length(b2sort))]
    b2med<- median(b2sort,na.rm=T)
    
    # b3temp<- exp(jstore)>1.5 
    b3temp<- exp(jstore)>1.2 
    b3temp<- colSums(b3temp,na.rm=T)/sum(where)
    b3sort<- sort(b3temp)
    b3top<- b3sort[ceiling(0.975*length(b3sort))]
    b3bot<- b3sort[ceiling(0.025*length(b3sort))]
    b3med<- median(b3sort,na.rm=T)
    
    # b4temp<- exp(jstore)>1.2
    # b4temp<- colSums(b4temp,na.rm=T)/sum(where)
    # b4sort<- sort(b4temp)
    # b4top<- b4sort[ceiling(0.975*length(b4sort))]
    # b4bot<- b4sort[ceiling(0.025*length(b4sort))]
    # b4med<- median(b4sort,na.rm=T)		
    
    # b1temp<- exp(jstore)<=0.4
    # b1temp<- colSums(b1temp,na.rm=T)/sum(where)
    # b1sort<- sort(b1temp)
    # b1top<- b1sort[ceiling(0.975*length(b1sort))]
    # b1bot<- b1sort[ceiling(0.025*length(b1sort))]
    # b1med<- median(b1sort,na.rm=T)
    
    # b2temp<- exp(jstore)<=0.8 & exp(jstore)>0.4
    # b2temp<- colSums(b2temp,na.rm=T)/sum(where)
    # b2sort<- sort(b2temp)
    # b2top<- b2sort[ceiling(0.975*length(b2sort))]
    # b2bot<- b2sort[ceiling(0.025*length(b2sort))]
    # b2med<- median(b2sort,na.rm=T)
    
    # b3temp<- exp(jstore)<=1.2 & exp(jstore)>0.8
    # b3temp<- colSums(b3temp,na.rm=T)/sum(where)
    # b3sort<- sort(b3temp)
    # b3top<- b3sort[ceiling(0.975*length(b3sort))]
    # b3bot<- b3sort[ceiling(0.025*length(b3sort))]
    # b3med<- median(b3sort,na.rm=T)
    
    # b4temp<- exp(jstore)>1.2
    # b4temp<- colSums(b4temp,na.rm=T)/sum(where)
    # b4sort<- sort(b4temp)
    # b4top<- b4sort[ceiling(0.975*length(b4sort))]
    # b4bot<- b4sort[ceiling(0.025*length(b4sort))]
    # b4med<- median(b4sort,na.rm=T)		
    
    
    ### Store results ###
    
    collapsed[y,]<- c(years[y],cmed,cbot,ctop)
    
    over[y,]<- c(years[y],omed,obot,otop)
    
    under[y,]<- c(years[y],umed,ubot,utop)
    
    jtemp<- apply(exp(jstore),2,median,na.rm=T) #calculate median of each column
    
    Meanitemp<- apply(exp(jstore),1,mean,na.rm=T) #store mean results for individual fisheries 

    Medianitemp<- apply(exp(jstore),1,median,na.rm=T) #store median results for individual fisheries 

    Sditemp<- apply(exp(jstore),1,sd,na.rm=T) #store median results for individual fisheries 
    
    isort<- t(apply(exp(jstore),1,sort))
    
    top<- isort[,ceiling(0.975*dim(isort)[2])]
    
    bot<- isort[,ceiling(0.025*dim(isort)[2])]
    
    
    
    ###Store results for individual fisheries###
    
    indbox=boxplot(t(jstore),plot=F)
    
    indtemp[where,1]<- as.character(core$id[where])
    indtemp[where,2]<- years[y]
    indtemp[where,3]<- Medianitemp
    indtemp[where,4]<- Meanitemp
    indtemp[where,5]<- top
    indtemp[where,6]<- bot
    indtemp[where,7]<- exp(indbox$stats[4,])
    indtemp[where,8]<- exp(indbox$stats[2,])
    indtemp[where,9]<- Sditemp
    
    
    jsort<- sort(jtemp) #Sort medians
    
    medians[y,]<- jsort
    
    bin[,1]<- c(b1med,b1bot,b1top)
    bin[,2]<- c(b2med,b2bot,b2top)
    bin[,3]<- c(b3med,b3bot,b3top)
    # bin[,4]<- c(b4med,b4bot,b4top)
    
    
    med<- mean(jsort,na.rm=T) #mean of medians
    top<- jsort[ceiling(0.975*length(jsort))] #upper CI
    bot<- jsort[ceiling(0.025*length(jsort))] #lower CI
    box<- boxplot(jsort,plot=F)
    
    
    
    c[y,]<- c(years[y],med,top,bot,box$stats[4],box$stats[2])
    # if (y==(length(years)-17))
    # {
    # show(dim(jstore))
    # quartz()
    # a=boxplot(t(jstore))
    # show(a)		
    # }
  }
  
  output<- list(Median=c,Individuals=indtemp,Collapsed=collapsed,Over=over,Under=under,Bin=bin,MedianDistribution=medians,DynamicDistribution=DistStorage)
  
  return(output)	
}


