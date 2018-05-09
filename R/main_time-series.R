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

# Auxiliary functions
sourceDirectory("./R/aux_functions", modifiedOnly=FALSE)


###
### Data
###

 trap_name <- c("GOM")

 # Creating output folder for each sediment trap
 if (!file.exists(paste("output/",trap_name,sep=""))){ dir.create(paste("output/",trap_name,"/",sep=""))}

 
  # Time-series data ****************************************************************** still define rules of next_open
  if(trap_name == "GOM"){
    raw_data <- get_gom_sst(overwrite = F)
    length(raw_data[,1])
    raw_data<- raw_data[86:nrow(raw_data),] # removing first part of GOM time-series with big gaps (next_open = 123 and 21)
    length(raw_data[,1])
  }# else{raw_data <- get_trap_data(overwrite = F)}
 
 
  # First differences time-series
  datafd <- get_first_diff(raw_data, trap_name, overwrite = F)
  
  # Data to use
  data <- raw_data[,-c(2:7)] # or <- datafd
  data$open <- ymd(data$open) # transforming to date format (lubridate pkg)

  # Plot scatterplot of everything
  # ggpairs(data[,2:ncol(data)])
  

###
### Analysis
###

  # Testing if time-series are stationary (see output/trap_name/stationary_test.csv)
  stationary <- test_stationary(data, trap_name, overwrite = F)
  
  # Cross-Correlation (lags)
  lags_correlation <- cross_correlation(data, trap_name, overwrite = F)
  
  # Splines and species seasonality
  splines <- seasonal_splines(DataSeries = data, DateFormat='%d/%m/%y', SavePlots=T)
  
  # Surrogates: randomization of residuals, summed with splines to generate null series (nreps = 500, i.e., 500 null series)
  corr_surrog <- surrogates_splines(data, splines, trap_name, nreps=500, overwrite=T)     
  # corr_surrog[1:50,1:10]
  
  # Plotting correlations with surrogate distribution (box-plots)
  plot_correlation_surrogates(corr_surrog, trap_name)
  


 