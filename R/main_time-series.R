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
library(rEDM)
library(zoo)

# Auxiliary functions
sourceDirectory("./R/aux_functions", modifiedOnly=FALSE)

# vignette("rEDM-tutorial", package="rEDM")
  
###### DATA ######

  # GOM_NAs.csv : Gulf of Mexico sediment trap data with continuous time-steps and "NA" for intervals without data, see "README.txt" in data folder
  data_na <- read.csv("data/GOM/GOM_NAs.csv", header = T, na = "NA")
  # changing format to date
  data_na$open <- dmy(data_na$open)
  data_na$close <- dmy(data_na$close)
  # transfoming columns into numeric
  data_na[, c(5:ncol(data_na))] <- sapply(data_na[, c(5:ncol(data_na))], function(x) as.numeric(as.character(x)))
  # warnings() OK! 'NA_resolution' and 'NA_gap' transformed into NA
  
  # Calculating the distribution of sums and differences between consecutives samples (less than 'days_closed' apart)
  distrib_list <- get_distrib_diffs(data_na, trap_name,  days_closed=7 ,overwrite = F)
  # Filling in the NA_resolutions
  data_use <- estimate_na_resolution(data_na, distrib_list, overwrite = F)
    
    
###### ANALYSIS ######
    
  data_ts <- ts(data_use[,4:19])
  # plot(data_ts[,c(2:11)])
  # plot(data_ts[,c(12:16)])
  
  # Embedding dimension: output/embedding_plots
  embed <- embed_dim(data_ts, emb_dim = 20, trap_name, overwrite = F) # uses simplex function rEDM, Rafa code
 
  # Simplex prediction and prediction decay: output/simplex_plots
  simplex_plot(data_ts, emax = embed$emax_eye, trap_name, overwrite = F)    
  
  # S-maps (theta): output/smap_plots
  # Red noise vs. nonlinear deterministic behaviour: if forecast skill increases for theta > 0, then the results are suggestive of nonlinear dynamics
  smap_plots(data_ts, emax = embed$emax_eye, trap_name, overwrite = F)
  
  
  ####################################################################################################
  # Convergent Cross Mapping (CCM)
  # Tutorial: https://mathbio.github.io/edmTutorials/ccm.html
  source("https://raw.githubusercontent.com/mathbio/edmTutorials/master/utilities/make_block.R")
  
  # max_lag is the optimal embedding dimension
  ssp_x <- c("G_ruber_pink")
  ssp_y <- c("G_ruber_white")
  
  embed[which(embed$species %in% c(ssp_x, ssp_y)),] 
  X <- c(data_ts[,which(colnames(data_ts) == ssp_x)])
  Y <- c(data_ts[,which(colnames(data_ts) == ssp_y)])
  paste("Total flux",ssp_x, sum(X, na.rm=T), "shells")
  paste("Total flux",ssp_y, sum(Y, na.rm=T), "shells")
  
  # Normalizing to mean 0 and unit variance 1
  X <- scale(X)
  Y <- scale(Y)
    
  fit<-lm(Y ~ X)
  plot(X,Y,main='Correlation (X,Y)')
  abline(0,1,lty=2)
  abline(fit$coefficients[1],fit$coefficients[2])
  legend(x = "bottomright", legend = paste('r =',round(cor(X,Y, use = "complete.obs"),2)))
  
  XY<-as.data.frame(cbind(X,Y)) 
  Shadow_MXY <- make_block(XY,max_lag = 2) # emax_eye = 2
  names(Shadow_MXY) <- c("time","X","X_1","Y","Y_1")
  head(Shadow_MXY)
  Shadow_MX<-Shadow_MXY[,2:3]
  Shadow_MY<-Shadow_MXY[,4:5]
  
  
  #####
  ##### X predicting Y (X cross-map Y)
  #####
  
  lib <- c(1, NROW(Shadow_MXY))  # cross-mapping process starts by finding (mapping) the simplex_Mx onto MY, creating the simplex_My. Note that simplex_My has the same indexes as simplex_Mx. 
  block_lnlp_output_XY <- block_lnlp(Shadow_MXY, lib = lib, pred = lib, columns = c("X","X_1"), target_column = "Y", stats_only = FALSE, first_column_time = TRUE)
  observed_all_Y <- block_lnlp_output_XY$model_output[[1]]$obs
  predicted_all_Y <- block_lnlp_output_XY$model_output[[1]]$pred
  pred_obs_Y<-as.data.frame(cbind(predicted_all_Y,observed_all_Y)) #
  colnames(pred_obs_Y)<-c('Predicted Y','Observed Y')
  head(pred_obs_Y)
  cor(observed_all_Y, predicted_all_Y,use = "complete.obs")
  
  # Plot
  fit_YX<-lm(predicted_all_Y ~ observed_all_Y)
  plot_range <- range(c(observed_all_Y, predicted_all_Y), na.rm = TRUE)
  plot(x=observed_all_Y,y=predicted_all_Y, xlim = plot_range, ylim = plot_range, xlab = "Observed Y",
       ylab = "Predicted Y", main='X cross-map Y')
  abline(fit_YX$coefficients[1],fit_YX$coefficients[2])
  abline(0,1,lty=2)
  legend(x = "topleft", legend = paste('r =',round(cor(observed_all_Y, predicted_all_Y,use = "complete.obs"),2)),inset = 0.02,col = 'black')
  
  #####
  ##### Y predicting X (Y cross-map X)
  #####
  
  lib<-lib <- c(1, NROW(Shadow_MXY))
  block_lnlp_output_YX <- block_lnlp(Shadow_MXY, lib = lib, pred = lib, columns = c("Y", "Y_1"),target_column = "X", stats_only = FALSE, first_column_time = TRUE)
  observed_all_X <- block_lnlp_output_YX$model_output[[1]]$obs
  predicted_all_X <- block_lnlp_output_YX$model_output[[1]]$pred
  pred_obs_X<-as.data.frame(cbind(predicted_all_X,observed_all_X))
  colnames(pred_obs_X)<-c('Predicted X','Observed X')
  head(pred_obs_X)
  cor(observed_all_Y, predicted_all_Y,use = "complete.obs")
  
  # Plot
  fit_XY<-lm(predicted_all_X ~ observed_all_X)
  plot_range <- range(c(observed_all_X, predicted_all_X), na.rm = TRUE)
  plot(observed_all_X, predicted_all_X, xlim = plot_range, ylim = plot_range,
       xlab = "Observed X", ylab = "Predicted X", main='Y cross-map X')
  abline(fit_XY$coefficients[1],fit_XY$coefficients[2])
  abline(0,1,lty=2)
  legend(x = "bottomright", legend = paste('r =',round(cor(observed_all_Y, predicted_all_Y,use = "complete.obs"),2)),inset = 0.02,col = 'black')
  
  #####
  #####   Convergent Cross Mapping (CCM)
  #####

    # cross map from X to Y
  X_xmap_Y<- ccm(XY, E = 2, lib_column = "X", target_column = "Y",
                 lib_sizes = seq(10, 170, by = 10), num_samples = 100, random_libs = TRUE,
                 replace = TRUE)
  # cross map from Y to X
  Y_xmap_X<- ccm(XY, E = 2, lib_column = "Y", target_column = "X",
                 lib_sizes = seq(10, 130, by = 10), num_samples = 100, random_libs = TRUE,
                 replace = TRUE)
  
  # mean values
  X_xmap_Y_means <- ccm_means(X_xmap_Y)
  Y_xmap_X_means <- ccm_means(Y_xmap_X)
  
  # plot graphs
  plot(X_xmap_Y_means$lib_size, pmax(0, X_xmap_Y_means$rho), type = "l", col = "red",
       main='Two Species', xlab = "Library Size (L)",
       ylab = "Cross Map Skill (Pearson rho)", ylim = c(0,1))
  lines(Y_xmap_X_means$lib_size, pmax(0, Y_xmap_X_means$rho), col = "blue")
  legend(x = "topleft", legend = c("X_xmap_Y", "Y_xmap_X"), col = c("red", "blue"),
         cex=1.1,lwd=2, inset = 0.02)
  