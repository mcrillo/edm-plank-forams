rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./time-series-forams/")

# Libraries
# source("R/library.R")
library(R.utils)   # sourceDirectory function
library(lubridate) # date calculations 
library(GGally)    # ggpairs function
library(tseries)   # tests series stationary
library(ggplot2)   # plots
library(reshape)   # function melt
library(corrplot)  # matrix correlation plot

# Auxiliary functions
sourceDirectory("./R/aux_functions", modifiedOnly=FALSE)


###
### Data
###

 trap_name <- c("GOM") # Gulf of Mexico sediment trap

 # Creating output folder for this sediment trap
 if (!file.exists(paste("output/",trap_name,sep=""))){ dir.create(paste("output/",trap_name,"/",sep=""))}

 
  # Time-series data 
  if(trap_name == "GOM"){
    raw_data <- get_gom_sst(overwrite = F) # gets temperature data for the GOM time-series
    raw_data<- raw_data[86:nrow(raw_data),] # removing first part of GOM time-series with big gaps (next_open = 123 and 21)
  }# Other sediment traps: else{raw_data <- get_trap_data(overwrite = F)} # still write "get_trap_data" function
 
 
  # Data to use
  data <- raw_data[,-c(2:7)] # or <- datafd
  data$open <- ymd(data$open) # transforming to date format (lubridate pkg)

  # Plot scatterplot of everything
  # ggpairs(data[,2:ncol(data)])
  
  # First differences time-series
  datafd <- get_first_diff(raw_data, trap_name, days_closed = 10,overwrite = F)
  # days_closed : maximum number of days inbetween samples,if that gap between two samples is bigger than 10 days (next_open > 10 days), then the difference between these two samples is not included

  
#########################
### Relative abundace ###
#########################
  
  # Relative abundance
  datarel <- get_relat_abund(raw_data, trap_name, overwrite=F)
  
  datarel_cor <- cormatrix(datarel[,-c(1:2)], cormethod = c("kendall"), conf.level = 0.95)
  datarel_mcor <- datarel_cor$corr
  colnames(datarel_mcor) <- rownames(datarel_mcor)<- colnames(datarel)[-c(1:2)]
  corrplot(datarel_mcor)
  
  # Plot total abundance through time
  p <- ggplot(datarel, aes(x=cum_days, y=total_abund)) + geom_line() + geom_point() +
    labs(y = "Flux (#shells * m-2 * day-1)", x = "Cumulative Days")  
  pdf(file =  paste("output/",trap_name,"/total_abund_",trap_name,".pdf",sep=""), width=12, height=8, paper = "special")
    print(p)
  dev.off()
  
  total_ssp_abund <- data.frame(total_ssp_abund = colSums(data[,-c(1,2)]))
  
  
############################
### Correlation analysis ###
############################

  # Testing if time-series are stationary (see output/trap_name/stationary_test.csv)
  stationary <- test_stationary(data, trap_name, overwrite = F)
  
  # Cross-Correlation (lags)
  lags_correlation <- cross_correlation(data, trap_name, overwrite = F)
  
  # Splines and species seasonality
  splines <- seasonal_splines(DataSeries = data, DateFormat='%d/%m/%y', SavePlots = F, overwrite = F)
  
  # Correlation and Surrogates
  corr_method <- c("kendall")
  corr_surrog <- corr_surrogates(data, splines, trap_name, nreps=500, corr_method, overwrite = T)   
  # Surrogates: randomization of residuals, summed with splines to generate null series (nreps = 500 null series, columns V1 - V500)
  # corr_surrog[1:50,1:15]
  
  # Plotting box-plots: correlations with surrogate distribution for each focal variable
  corr_surrogates_boxplots(corr_surrog, trap_name, overwrite=T)
   

#####################################
### EDM analysis: GOM time-series ###
#####################################
  
  library(rEDM)
  library(zoo)
  
###### DATA ######

  # GOM_NAs.csv : Gulf of Mexico sediment trap data with continuous time-steps and "NA" for intervals without data, see "README.txt" in data folder
  data_na <- read.csv("data/GOM/GOM_NAs.csv", header = T, na = "NA")
  # changing format to date
  data_na$open <- dmy(data_na$open)
  data_na$close <- dmy(data_na$close)
  # tranfoming columns into numeric, NA_resolution and NA_gap will be transformed into NA (warnings())
  data_na[, c(5:ncol(data_na))] <- sapply(data_na[, c(5:ncol(data_na))], function(x) as.numeric(as.character(x)))
  
  # Calculating the distribution of sums and differences between consecutives samples (less than 'days_closed' apart)
  distrib_list <- get_distrib_diffs(data_na, trap_name,  days_closed=7 ,overwrite = F)
  # Filling in the NA_resolutions
  data_use <- estimate_na_resolution(data_na, distrib_list, overwrite = F)
    
    
###### ANALYSIS ######
    
  data_ts <- ts(data_use[,4:19])
  # plot(data_ts[,c(2:11)])
  # plot(data_ts[,c(12:16)])
  
  # Embedding dimension plot (simplex function rEDM) # Rafa
  embed_max <- embed_plots(data_ts, emb_dim = 20, trap_name, overwrite = F)   
  names(embed_max) <- c("species","emax_auto")
  # Optimizing embeding dimension by eye and hand:
  embed_max <- cbind(embed_max, emax = c(3,3,5,2,2,2,1,2,2,5,3,10,3,5,5))
  # worse eye-optimization: G_ruber_pink, N_dutertrei, G_truncatulinoides
  
  # Simplex prediction and prediction decay
  simplex_plot(data_ts, emax = embed_max$emax, trap_name, overwrite = F)    
  
  # S-maps: red noise vs. nonlinear deterministic behaviour
  # If forecast skill increases for theta > 0, then the results are suggestive of nonlinear dynamics
  smap_plots(data_ts, emax = embed_max$emax, trap_name, overwrite = F)
  
  
  
    