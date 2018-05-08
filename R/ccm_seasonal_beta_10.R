ccm_seasonal_beta_10 <- function(FileFolder=getwd(),
                                 FileName='plankton_final.csv',
                                 TargetColumn='G_bul',
                                 DateColumn='time',
                                 num_surr=500,
                                 T_period=52,
                                 seed=12345,
                                 significance=0.05,
                                 MaxE=10,
                                 XFilter=0,
                                 image.width=40,
                                 image.height=20,
                                 Seasonality=T,
                                 PrintSeason=T,
                                 DateFormat=F,
                                 UseAIC=F,
                                 PreAnalysis=F,
                                 PlotEPS=F,
                                 TimeForPred=0){
  #FileFolder: patch of data file
  #FileName: dara file name with extension [*.csv]
  #TargetColumn: name of target variable(effect) [ex. 'vivax.cases']
  #DateColumn: name of date column [usually 'date']
  #num_surr: number of surrogates to build null hypothesis distribution
  #T_period: period of seasonality in units of sampling period  [ex. date in months and season of 12 months,T_period<-12]
  #seed: set seed for reproducible results
  #significance: significance level for hypothesis testing [usually significance<-0.05]
  #assumes all other columns other than 'TargetColumn' and 'DateColumn' are possible causes
  # MaxE: maximum embedding dimension to optimize attractor reconstruction
  # XFilter: number of times the moving average is applied
  # Seasonality: (TRUE) if there is seasonality in the data
  # PrintSeason: (TRUE) if you want to generate seasonality figures
  # DateFormat: (TRUE) if date is in standard form
  # UseAIC: (TRUE) if AIC criteria is neede for embedding optimization
  # PreAnalysis: (TRUE) if nonlinearity and convergence needs to be tested
  # PlotEPS: (TRUE) figures in EPS format file
  # TimeForPred: delay between cause and effect (use negative values)
  #*********************************
  # function outputs
  #- curves comparing the real data with one surrogate for each variable
  #- boxplot with surrogate distribution from all possible cause variable.
  #ccm value for real data represented by an empty red circle (p>significance) and
  #full red circle if (p<=significance)
  #- table containing all ccm values for surrogate and real data (first line)
  #- table showing the optimal embeding and p-values for each variable
  #- a date/time stamp is given for each generated output, therefore it will not
  #overwrite previous runs
#*********************************
  #this function uses libraries
  #library('rEDM')
  #library('plotly')
  #library('tidyr')
  #library('dplyr')
  #library('ggplot2')
  #library('reshape2')
#begin function*********************************
  set.seed(seed)
  #_________________________________________________________________________________
  # setting files parameters
  dir.create(paste(FileFolder,'/results_pre-processing/',sep=""),showWarnings = FALSE)
  FileNameSave<-substr(FileName,1,nchar(FileName)-4)
  setwd(FileFolder)
  DatePrint<-paste(as.character(format(Sys.time(), "%Y")),as.character(format(Sys.time(), "%b")),as.character(format(Sys.time(), "%d")),sep="-")
  TimePrint<-gsub(":", "-", as.character(format(Sys.time(), "%X")))
  TimeStamp<-paste(DatePrint,TimePrint,sep="__")
  NewFolder<-paste(FileFolder,'/results/',TimeStamp,'__',FileNameSave,' lag=',TimeForPred,' target=',TargetColumn,sep="")
  dir.create(NewFolder)
  #_________________________________________________________________________________
  #FileNameSave<-substr(FileName,1,nchar(FileName)-4)
  ##time_series: all data with appropriate headings
  #setwd(FileFolder)
  #timestamp<-as.numeric(Sys.time())
  #dir.create(paste(FileFolder,'/results/',sep=""),showWarnings = FALSE)
  #NewFolder<-paste(FileFolder,'/results/',timestamp,'__',FileNameSave,sep="")
  #dir.create(NewFolder)


  ##pre-processing steps
  all.data <- read.csv(FileName,header = TRUE)
  variables<-names(all.data)
  #distinguishing variables
  #v.target<-variables[which(variables==TargetColumn)]
  v.target.index<-which(variables==TargetColumn)
  v.date.index<-which(variables==DateColumn)
  v.drivers<-variables[-c(v.target.index,v.date.index)]


  VPreAnalysis<-variables[-v.date.index]
  VPreAnalysisIndex<-which(variables %in% VPreAnalysis)


  v.drivers.index<-which(variables %in% v.drivers)
  if (DateFormat==T){
  v.date.data<-as.Date(strptime(all.data[,v.date.index],format="%d/%m/%Y"))
  }
  else{
    v.date.data<-all.data[,v.date.index]
  }

  #matrices for results presentation
  final_ccm_surr<-matrix(NA,num_surr,length(v.drivers))
  final_ccm_original<-matrix(NA,1,length(v.drivers))
  final_boxplot_result<-matrix(NA,num_surr+1,length(v.drivers))
  E_star <- matrix(NA,1,length(v.drivers))
  p.value<-matrix(NA,1,length(v.drivers))
  significative<-matrix(NA,length(v.drivers),2)
  #progress bar
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = length(v.drivers), width = 300)
#*************************************************
  # pre analysis test
  if (PreAnalysis==T){
    FolderPreAnalysis<-paste(NewFolder,'/PreAnalysis/',sep="")
    dir.create(FolderPreAnalysis)
    for (i in 1:length(VPreAnalysis)){
      VAuxData<-all.data[,VPreAnalysisIndex[i]]
      SimplexOutput <- simplex(VAuxData,E=1:MaxE)
      OptimalE<-which.max(SimplexOutput$rho[2:length(SimplexOutput)])+1
      TimeToPrediction <- simplex(VAuxData, E = OptimalE, tp = 1:10)
      NonLinearity <- s_map(VAuxData, E = OptimalE)

      # plot simplex
      png(paste(FolderPreAnalysis,'/_simplex_',VPreAnalysis[i],'.png',sep=''),width=800,height=600)
      plot(SimplexOutput$E,SimplexOutput$rho,type='l',xlab = "Embedding Dimension (E)",ylab = "Forecast Skill (rho)",ylim=c(0, 1),main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      lines(OptimalE,SimplexOutput$rho[OptimalE],type='p',col='red',pch = 19)
      dev.off()
      # plot PredictionToPrediction
      png(paste(FolderPreAnalysis,'/_prediction_',VPreAnalysis[i],'.png',sep=''),width=800,height=600)
      plot(TimeToPrediction$tp,TimeToPrediction$rho,type='l',xlab = "Time to Prediction (tp)",ylab = "Forecast Skill (rho)",ylim=c(0, 1),main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      dev.off()
      # plot NonLinearity
      png(paste(FolderPreAnalysis,'/_nonlinearity_',VPreAnalysis[i],'.png',sep=''),width=800,height=600)
      plot(NonLinearity$theta,NonLinearity$rho,type='l',xlab = "Nonlinearity (theta)",ylab = "Forecast Skill (rho)",main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      dev.off()

      if (PlotEPS==T){
        NewFolderEPS<-paste(FolderPreAnalysis,'/eps',sep='')
        dir.create(NewFolderEPS)
        # plot simplex
        setEPS()
        postscript(paste(NewFolderEPS,'/','/_simplex_',VPreAnalysis[i],'.eps',sep=''),width=8.5,height=6)
        plot(SimplexOutput$E,SimplexOutput$rho,type='l',xlab = "Embedding Dimension (E)",ylab = "Forecast Skill (rho)",ylim=c(0, 1),main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        lines(OptimalE,SimplexOutput$rho[OptimalE],type='p',col='red',pch = 19)
        dev.off()
        # plot PredictionToPrediction
        setEPS()
        postscript(paste(NewFolderEPS,'/','/_prediction_',VPreAnalysis[i],'.eps',sep=''),width=8.5,height=6)
        plot(TimeToPrediction$tp,TimeToPrediction$rho,type='l',xlab = "Time to Prediction (tp)",ylab = "Forecast Skill (rho)",ylim=c(0, 1),main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        dev.off()
        # plot NonLinearity
        setEPS()
        postscript(paste(NewFolderEPS,'/','/_nonlinearity_',VPreAnalysis[i],'.eps',sep=''),width=8.5,height=6)
        plot(NonLinearity$theta,NonLinearity$rho,type='l',xlab = "Nonlinearity (theta)",ylab = "Forecast Skill (rho)",main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        dev.off()
      }

# **************************
# ARRUMAR TAMANHO DOS LABELS E NUMEROS DOS GRAFICOS
# IMPLEMENTAR ESTUDO CONVERGENCIA
# **************************

      #MaxLib<-floor(NROW(all.data)*0.9)
      #LibPre<-
      remove(VAuxData)
    }
  }
#*************************************************


#loop for causes variables (drivers)
for (i in 1:length(v.drivers)){
  #print(i)
  E_star_surr<-matrix(NA,num_surr,1)
    v.aux<-v.drivers[i]
    v.aux.index<-v.drivers.index[i]
    v.aux.data<-all.data[,v.aux.index]
    v.target.data<-all.data[,v.target.index]
    lib_ccm <- c(1,NROW(all.data))
#***************************************************



#***************************************************
    ##original data
    #optimal embeding
    block_temp <- select(all.data,v.target.index,v.aux.index)

    #Begin filtering - original data*************************************
    block_temp_aux<-block_temp

    block_temp[,1]<-FilterAverage(block_temp[,1],XFilter)
    block_temp[,2]<-FilterAverage(block_temp[,2],XFilter)
    #*
    #plot filtered original data
    #plot(block_temp_aux[,1])
    #plot(block_temp[,1],type='l')
    #readline(prompt="Press [enter] to continue")
    #*
    v.target.data_filtered<-FilterAverage(v.target.data,XFilter)
    #End filtering - original data*************************************
    out.temp <- do.call(
      rbind,
      lapply(1:MaxE, function(E_i){
        pred_ccm <- make_pred_nozero(v.target.data_filtered,E_i)
        #ccm
        ccm(block=block_temp,
            E=E_i,
            lib=lib_ccm,
            pred=pred_ccm,
            lib_sizes = NROW(block_temp),
            exclusion_radius=0,
            random_libs = FALSE,
            num_sample=1,
            tp = TimeForPred-1,
            lib_column = 1,
            target_column = 2)
      })
    )
    #return(out.temp)
    #*
    #(2)#plot - rho(E) - original data
    #return(out.temp)
    #par(mfrow=c(2,1))
    #plot(out.temp$rho,type='l')
    #plot(out.temp$rmse,type='l')

    #readline(prompt="Press [enter] to continue")

    #readline(prompt="Press [enter] to continue")
    #*
    #return(out.temp)


    #*AIC
    if (UseAIC==T){
    Vector_AUX<-matrix(NA,nrow(out.temp),2)
    Vector_AUX[,1]<-out.temp$rmse
    Vector_AUX[,2]<-out.temp$E
    out.temp$rho<-AIC_CCM(Vector_AUX[,1],Vector_AUX[,2],nrow(all.data)-Vector_AUX[,2])
    AuxPeak<-((out.temp$rho-min(out.temp$rho))-2)
    E_star[i]<-out.temp$E[min(which(AuxPeak<0))]
    }
    else{
      E_star[i] <- out.temp$E[which.max(out.temp$rho[2:MaxE])+1]
    }
    #*AIC

    #(2)#plot - rho(E) - original data
    SaveName<-paste(NewFolder,'/',FileNameSave,v.aux,'_Rho_E','.png',sep="")
    png(SaveName)
    plot(out.temp$rho,type='l',main=v.aux,xlab="E",ylab="rho")
    #plot(out.temp$rmse,type='l',main=v.aux,xlab="E",ylab="rho")


    dev.off()
    #readline(prompt="Press [enter] to continue")
    #*

    #*AIC
    #AuxPeak<-((out.temp$rho-min(out.temp$rho))-2)
    #E_star[i]<-out.temp$E[min(which(AuxPeak<0))]
    #return(E_star[i])
    #*AIC

    #E_star[i] <- out.temp$E[which.max(out.temp$rho)]


    pred_ccm <- make_pred_nozero(v.target.data,E_star[i])
    #ccm
    df.out.ccm <- ccm(block=block_temp,
                      E=E_star[i],
                      lib=lib_ccm,
                      pred = pred_ccm,
                      lib_sizes = NROW(block_temp),
                      exclusion_radius=0,
                      random_libs = FALSE,
                      num_sample=1,
                      tp = TimeForPred)
#***************************************************
#***************************************************
    ##surrogate data
    v.aux.surr <- make_surrogate_data(v.aux.data,
                                      method = "seasonal",
                                      T_period = T_period,
                                      num_surr = num_surr)
    #***********************
     #for (k in 1:dim(v.aux.surr)[2])
      # plot(v.aux.surr[1:100,k],type='l')
      #{
       #size.new<-dim(v.aux.surr)[1]-2
       #denoise<-matrix(NA,size.new,dim(v.aux.surr)[2])
       #series.aux<-v.aux.surr[,k]
       #at.aux<-size.new
      #for (h in 2:at.aux)
      #{
       # start<-h-1
        #finish<-h+1
        #v.aux.surr[h,k]<-mean(series.aux[start:finish])
      #}
      #plot(v.aux.surr[1:100,k],type='l')
     #}
    #***********************
    if (PrintSeason==T){
    #plot seasonality

      SaveName<-paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_seasonality_',v.aux,'.png',sep="")
      png(SaveName)
      plot(v.date.data,v.aux.data,type='l',xlab="date",ylab=v.aux,width=800,height=600)
      lines(v.date.data,rowMeans(v.aux.surr),col='red')
      dev.off()
    #ggplot(all.data, aes(v.date.data)) +
    #  geom_line(aes(y = v.aux.data, colour = 'original'),size=1)+
    #  geom_line(aes(y = rowMeans(v.aux.surr), colour = "seasonality"),size=1)+
    #  labs(x='date',y=v.aux)
    #ggsave(paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_seasonality_',v.aux,'.pdf',sep=""),
    #       width = image.width, height = image.height, units = "cm")
    }
    #optimal embeding
    df.out.ccm.surr <- data.frame()
    for (i_surr in 1:dim(v.aux.surr)[2]){
      remove(block_temp)
      block_temp <- data.frame(target.block=v.target.data,v.block=v.aux.surr[,i_surr])
      #Begin filtering - surrogate data*************************************
      block_temp[,1]<-FilterAverage(block_temp[,1],XFilter)
      block_temp[,2]<-FilterAverage(block_temp[,2],XFilter)
      #*
      #plot filtered surrogate data
      #plot(block_temp_aux[,1])
      #plot(block_temp[,1],type='l')
      #readline(prompt="Press [enter] to continue")
      #*
      #End filtering - surrogate data*************************************
      out.temp <- do.call(
        rbind,
        lapply(1:MaxE, function(E_i){
          pred_ccm <- make_pred_nozero(v.target.data_filtered,E_i)
          ccm(block=block_temp,
              E=E_i,
              lib=lib_ccm,
              pred=pred_ccm,
              lib_sizes = NROW(block_temp),
              exclusion_radius=0,
              random_libs = FALSE,
              num_sample=1,
              tp = TimeForPred-1,
              lib_column = 1,
              target_column = 2)
        })
      )

      #*AIC
      if (UseAIC==T){
        Vector_AUX<-matrix(NA,nrow(out.temp),2)
        Vector_AUX[,1]<-out.temp$rmse
        Vector_AUX[,2]<-out.temp$E
        out.temp$rho<-AIC_CCM(Vector_AUX[,1],Vector_AUX[,2],nrow(all.data)-Vector_AUX[,2])
        AuxPeak<-((out.temp$rho-min(out.temp$rho))-2)
        E_star_surr[i_surr]<-out.temp$E[min(which(AuxPeak<0))]
      }
      else{
        E_star_surr[i_surr] <- out.temp$E[which.max(out.temp$rho[2:MaxE])+1]
      }
      #*AIC

      #*AIC surrogate
      #Vector_AUX<-matrix(NA,nrow(out.temp),2)
      #Vector_AUX[,1]<-out.temp$rmse
      #Vector_AUX[,2]<-out.temp$E
      #out.temp$rho<-AIC_CCM(Vector_AUX[,1],Vector_AUX[,2],nrow(all.data)-Vector_AUX[,2])
      #*AIC surrogate

      #*AIC
      #AuxPeak<-((out.temp$rho-min(out.temp$rho))-2)
      #E_star_surr[i_surr]<-out.temp$E[min(which(AuxPeak<0))]
      #return(E_star[i])
      #*AIC






      #*
      #(3)#plot - rho(E) - surrogate data
      #plot(out.temp$rho,type='l')
      #readline(prompt="Press [enter] to continue")
      #*


      # E_star_surr[i_surr] <- out.temp$E[which.max(out.temp$rho)]
#***************************************************
##**## distribuicao E_star_surr
#***************************************************
      pred_ccm <- make_pred_nozero(v.target.data,E_star_surr[i_surr])
      out.temp <- ccm(block=block_temp,
                      E=E_star_surr[i_surr],
                      lib=lib_ccm,
                      pred = pred_ccm,
                      lib_sizes = NROW(block_temp),
                      exclusion_radius=0,
                      random_libs = FALSE,
                      num_sample=1,
                      tp = TimeForPred)
      df.out.ccm.surr <- df.out.ccm.surr %>% bind_rows(out.temp)
    } # end for i_surr
    #plot and save surrogate
#    ggplot(all.data, aes(v.date.data)) +
#      geom_line(aes(y = v.aux.data, color = 'original'),size=1)+
#      geom_line(aes(y = v.aux.surr[,1], colour = "surrogate"),size=1)+
#      labs(x='date',y=v.aux)
#    ggsave(paste(timestamp,v.aux,'_surrogate_02.png'))
   final_ccm_original[i]<-df.out.ccm$rho
   final_boxplot_result[2:nrow(final_boxplot_result),i]<-df.out.ccm.surr[,9]
   E_star_surr<-as.data.frame(E_star_surr)
   #*
   #*plot E* histogram for surrogate
   #ggplot(data=E_star_surr, aes(E_star_surr$V1)) +
   #geom_histogram(col="red",
   #                  fill="green",
   #                  alpha = .2) +
   #labs(title="Histogram for E Surrogates") +
   #labs(x="E", y="Count") +
   #xlim(c(min(E_star_surr$V1),max(E_star_surr$V1))) +
   #ylim(c(0,num_surr/3))
   #ggsave(paste(NewFolder,'/',timestamp,'_',FileNameSave,'_HistESurr_',v.aux,'.pdf',sep=''),
   #        width = image.width, height = image.height, units = "cm")
   #*
   #progress bar
   Sys.sleep(0.1)
   setWinProgressBar(pb, i, title=paste( round(i/length(v.drivers)*100, 0),"% complete"))
}#end loop variables
  close(pb)
  final_boxplot_result[1,]<-final_ccm_original
  colnames(final_boxplot_result)<-v.drivers
  #*csv with all values of rho for surrogates, first row is rho for the original data.
  write.csv(final_boxplot_result,paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_rho.csv',sep=''))
  #*
  #p.values
  for (i in 1:length(v.drivers.index)){
    p.value[i]<-(sum(final_boxplot_result[,i] > final_boxplot_result[1,i]) + 1) / (nrow(final_boxplot_result))
  }
  significative_logical<-p.value<significance
  total_sig<-sum(significative_logical, na.rm=TRUE)
  v.index.sig<-which(significative_logical==TRUE)
  rho.index.sig<-final_ccm_original[v.index.sig]
  if (total_sig>0){
  for (i in 1:total_sig){
  significative[i,1]<-v.index.sig[i]
  significative[i,2]<-rho.index.sig[i]
  }}
  results_final<-rbind(E_star,p.value)
  significative<-as.data.frame(significative)
  colnames(results_final)<-v.drivers
  rownames(results_final)<-c('E*','p-value')
  #*csv with p-values and E*
  write.csv(results_final,paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_results.csv',sep=''))
  #*
  highlight <- c(TRUE, rep(FALSE, nrow(final_boxplot_result) - 1L))  # tag first row as interesting
  highlight<-rep(highlight,5)
  df.2 <- melt(final_boxplot_result)  # convert df to long format
  colnames(df.2)<-c('n','variable','rho')

  #png(paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_boxplot.png',sep=''))
  ggplot(subset(df.2, !highlight), aes(x=variable, y=rho)) +
    geom_boxplot()  +
    geom_point(data=subset(df.2, highlight),
               aes(x=variable, y=rho),
               color="red", size=4, shape=21, stroke = 1, fill=NA) + {if (total_sig>0) geom_point(data=significative,
                                                                   aes(x=significative$V1, y=significative$V2),
                                                                   color = 'red', size=4, shape=21, fill='red') }
  #dev.off()
  #saving boxplot figure


  #plot(SimplexOutput$E,SimplexOutput$rho,type='l',xlab = "Embedding Dimension (E)",ylab = "Forecast Skill (rho)",ylim=c(0, 1),main=VPreAnalysis[i],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #lines(OptimalE,SimplexOutput$rho[OptimalE],type='p',col='red',pch = 19)

  ggsave(paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_boxplot.png',sep=''),
         width = image.width, height = image.height, units = "cm")
  #save file with parameters e main results
  fileConn<-file(paste(NewFolder,'/',TimeStamp,'_',FileNameSave,'_parameters.txt',sep=""))
  writeLines(c(paste('Date and Time:', Sys.time()),
               paste('File:',FileName),
               paste('TargetColumn:', TargetColumn),
               paste('DateColumn:', DateColumn),
               paste('num_surr:', num_surr),
               paste('T_period:', T_period),
               paste('seed:', seed),
               paste('significance:', significance),
               paste('MaxE:', MaxE),
               paste('XFilter:', XFilter),
               paste('image.width:', image.width),
               paste('image.height:', image.height)), fileConn)
  close(fileConn)
}