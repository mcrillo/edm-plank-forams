
seasonal_splines <- function(DataSeries,DateFormat, SavePlots, overwrite){
  #---
  # Author: Brenno Cabella (SeasonalSplines.R modified by Marina Rillo)
  # This code extracts the sazonality from time series using splines.
  # Inputs:   DataSeries - Time-series with only one date ("open") column
  #           DateFormat - Format of date (usually '%d/%m/%Y')
  #           SavePlots - Save sazonality results as figures (T or F)
  # Outputs:  Figure if SavePlots=T
  #           csv file with the seasonality for all variables in DataFile   
  #---
  
if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/splines_GOM.csv",sep=''))){
    
  # read data
  DateCol=1 # "open" column
  # names variables
  VariablesNames<-names(DataSeries)
  # attribute date format to date column
  DataSeries[,DateCol]<-as.Date(DataSeries[,DateCol], DateFormat)
  # save original date format for plot
  OriginalDate<-DataSeries[,DateCol]
  # change date format to track days and months only
  AuxDate<-format(DataSeries[,DateCol],"%m/%d")
  DataSeries[,DateCol]<-format(as.Date(DataSeries[,DateCol], DateFormat), format='%m/%d')
  # create output matrix
  Seasonal<-matrix(NA,nrow(DataSeries),ncol(DataSeries))
  # agregate same month-day, taking the mean of species' flux in the same month-day over the years
  AgDados<-aggregate(DataSeries[,-DateCol], by=list(DataSeries[,DateCol]), FUN=mean, na.rm=TRUE)
  # weight for the splines considering the repetitions of each month-day
  WDados <- table(DataSeries[,DateCol])
  # match the dates in the original data with the aggregated ones
  MatchDate<-match(AuxDate,AgDados[,1])
  
  # loop for variables
  for (i in 2:ncol(DataSeries)){
    # replicate 3 times weights and smooth to avoid discontinuities
    # 3 times so that the middle series (the one we are using) does not have weird discontinuities in the beginning and the end
    AgDadosRep<-rep(AgDados[,i],3)
    spline <- smooth.spline(AgDadosRep,w=rep(WDados,3))
    
    # str(spline)
    # all(spline$data$y == AgDadosRep) # spline$data is the original data
    SeasonalDataRep<-spline$y # estimated seasonality
    
    # take the middle part of the spline result (all 3 parts are identical)
    SeasonalDataAux<-SeasonalDataRep[(floor(length(SeasonalDataRep)/3)+1):floor((2*length(SeasonalDataRep)/3))]
    # store result in output matrix
    Seasonal[,i]<-SeasonalDataAux[MatchDate]
    
    # test of how seasonal a species is: correlation of raw data with spline
    # SspSeasonal <- rbind(SspSeasonal, data.frame(variable = VariablesNames[i], correlation = cor(spline$data$y, spline$y)))
    # write.csv(SspSeasonal,"output/GOM_splines/seasonal_species_corr.csv",row.names = F)
    # cor(DataSeries[,i],SeasonalDataAux[MatchDate])
 
    # save plot
    if (SavePlots==T){
      if (!file.exists(paste("output/",trap_name,"/splines_plots",sep=""))){ dir.create(paste("output/",trap_name,"/splines_plots/",sep=""))}
      setEPS()
      postscript(paste("output/",trap_name,"/splines_plots/",VariablesNames[i],'_GOM.eps',sep=''),width=8.5,height=6)
      plot(OriginalDate,DataSeries[,i],ylab=VariablesNames[i],xlab='date')
      lines(OriginalDate,Seasonal[,i],col='blue',lwd=2)
      dev.off()
    } # if
  } # for
  
  # store original date in output
  Seasonal[,1]<-format(as.Date(OriginalDate,'%y-%m-%d'),'%Y-%m-%d')
  # insert columns names
  colnames(Seasonal)<-VariablesNames
  # write csv output file
  
  write.csv(Seasonal,paste("output/",trap_name,"/splines_GOM.csv",sep=''),row.names = F)
  
  splines <- as.data.frame(Seasonal, stringsAsFactors = FALSE)
  splines[,2:ncol(splines)]<- sapply(splines[,2:ncol(splines)], as.numeric)
  splines$open <- ymd(splines$open) # transforming to date format
  # data.frame(names(gom_columns), names(splines))
  
  
  # Seasonality test: correlation between a time-series and its spline
  seasonality_species <- data.frame(species = rep(NA,ncol(splines)), corr_spline = rep(NA,ncol(splines)))
  for (i in 2:ncol(splines)){
    seasonality_species[i,"species"] <- names(splines)[i]
    seasonality_species[i,"corr_spline"] <- cor(splines[,i], data[,i])
  }
  seasonality_species <- seasonality_species[-1,]
  write.csv(seasonality_species, paste("output/",trap_name,"/seasonality_GOM.csv",sep=''),row.names = F)
  
  
  return(splines)
  
}else{
  
  splines <- read.csv(paste("output/",trap_name,"/splines_GOM.csv",sep=''), header = TRUE)
  return (splines)
}
  
}


