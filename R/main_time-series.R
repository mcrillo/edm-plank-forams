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

 trap_name <- c("GOM") # Gulf od Mexico sediment trap

 # Creating output folder for this sediment trap
 if (!file.exists(paste("output/",trap_name,sep=""))){ dir.create(paste("output/",trap_name,"/",sep=""))}

 
  # Time-series data 
  if(trap_name == "GOM"){
    raw_data <- get_gom_sst(overwrite = F) # gets temperature data for the GOM time-series
    raw_data<- raw_data[86:nrow(raw_data),] # removing first part of GOM time-series with big gaps (next_open = 123 and 21)
  }# Other sediment traps: else{raw_data <- get_trap_data(overwrite = F)} # still write "get_trap_data" function
 
 
  # First differences time-series
  datafd <- get_first_diff(raw_data, trap_name, overwrite = F)
  
  # Data to use
  data <- raw_data[,-c(2:7)] # or <- datafd
  data$open <- ymd(data$open) # transforming to date format (lubridate pkg)

  # Plot scatterplot of everything
  # ggpairs(data[,2:ncol(data)])
  

############################
### Correlation analysis ###
############################

  # Testing if time-series are stationary (see output/trap_name/stationary_test.csv)
  stationary <- test_stationary(data, trap_name, overwrite = F)
  
  # Cross-Correlation (lags)
  lags_correlation <- cross_correlation(data, trap_name, overwrite = F)
  
  # Splines and species seasonality
  splines <- seasonal_splines(DataSeries = data, DateFormat='%d/%m/%y', SavePlots = T, overwrite = F)
  
  
  # Correlation and Surrogates
  corr_method <- c("kendall")
  corr_surrog <- corr_surrogates(data, splines, trap_name, nreps=500, corr_method, overwrite = T)   
  # Surrogates: randomization of residuals, summed with splines to generate null series (nreps = 500 null series, columns V1 - V500)
  # corr_surrog[1:50,1:15]
  
  # Plotting box-plots: correlations with surrogate distribution for each focal variable
  corr_surrogates_boxplots(corr_surrog, trap_name, overwrite=T)
   

####################
### EDM analysis ###
####################
  library(rEDM)
  
  # Gulf of Mexico sediment trap with continuous time-steps and "NA" for intervals without data
  data_na <- read.csv("data/GOM/GOM_NA_gap.csv", header = T, na = "NA")
  data_ts <- ts(data_na[,4:19])
  # plot(data_ts[,c(2:11)])
  # plot(data_ts[,c(12:16)])
  
  # Embedding dimension plot (simplex function rEDM) # Rafa
  embed_max <- embed_plots(data_ts, emb_dim = 20, trap_name, overwrite =F)   
  
  # Simplex prediction and prediction decay
  simplex_plot(data_ts, embed_max, trap_name, overwrite = F)    
  
  # S-maps: red noise vs. nonlinear deterministic behaviour
  # If forecast skill increases for θ>0, then the results are suggestive of nonlinear dynamics
  for (i in 2:ncol(data_ts)){
  X <- as.data.frame(data_ts)[,i]

  smap_output <- s_map(X, lib = c(1, NROW(X)), pred = c(1, NROW(X)), E = embed_max[i-1,"emax"], silent = T)
  
  savename<-paste("output/",trap_name,"/smap_plots/theta_",colnames(data_ts)[i],".png",sep="")
  png(savename, width = 800, height = 600)
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
  plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", 
       ylab = "Forecast Skill (rho)")
  dev.off()
  }