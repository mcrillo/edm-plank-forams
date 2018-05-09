rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./time-series-forams/")

# Libraries
# source("R/library.R")
library(R.utils)   # sourceDirectory function
library(lubridate) # date calculations 
library(GGally)    # ggpairs function
library(tseries)   # tests series stationary

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
  
  # Splines
  splines <- seasonal_splines(DataSeries = data, DateFormat='%d/%m/%y', SavePlots=F)
  
  # Seasonality test: correlation between a time-series and its spline
  seasonality_species <- data.frame(species = rep(NA,ncol(splines)), corr_spline = rep(NA,ncol(splines)))
  for (i in 2:ncol(splines)){
    seasonality_species[i,"species"] <- names(splines)[i]
    seasonality_species[i,"corr_spline"] <- cor(splines[,i], data[,i])
  }


  # Surrogates 

  # Randomization of residuals of (data - splines)
  # data.frame(data = gom_columns[,i], spline = splines[,i], resid = gom_columns[,i] - splines[,i])
  nreps = 500
  random_resid <- list()
  random_resid <- lapply(2:ncol(splines), function(i) replicate(nreps, sample(c(gom_columns[,i] - splines[,i]), size = length(splines[,i]), replace = FALSE, prob = NULL)))
  names(random_resid) <- names(splines)[2:ncol(splines)]
  str(random_resid)
  # var(random_resid[[i]][,50]) # double check if variance is always the same, i.e., if vectors are all identical just with values in diferent order (randomized)

  # Generating surrogates for each species: adding random residuals to splines (i.e., seasonal time-series)
  # data.frame(resid = random_resid[[i]][,50], spline = splines[,i],  surrogate = random_resid[[i]][,50]+ splines[,i])
  surrogates <- lapply(1:length(random_resid), function(i){random_resid[[i]] + splines[,i+1]}) # i+1 because first column of splines is "open"
  names(surrogates) <- names(random_resid)
  str(surrogates)
  # all(round(surrogates[[i]][,j] - random_resid[[i]][,j] - splines[i+1],4)==0) # double check: this must be all 0 for any value of i in [1,15] and j[1,500]

  # Calculate correlations column by columns of each variable pair
  corr_pairs <- data.frame()
  for(i in 1:length(comb_two[1,])){
  
    col_var1 <- which(names(gom_columns)==names(surrogates)[(comb_two[1,i])])
    col_var2 <- which(names(gom_columns)==names(surrogates)[(comb_two[2,i])])
  
    # correlation of each variable with its own spline (measure of seasonality of the time-series)
    season_var1 <- cor(gom_columns[,col_var1], splines[,col_var1])
    season_var2 <-  cor(gom_columns[,col_var2], splines[,col_var2])
    
    # correlations
    corr_data <- cor(gom_columns[,col_var1], gom_columns[,col_var2]) # red line in the hist
    corr_splines <- cor(splines[,col_var1], splines[,col_var2])  # black line in the hist
    corr_resid <- cor( (gom_columns[,col_var1]-splines[,col_var1]), (gom_columns[,col_var2]-splines[,col_var2]) )  # blue line in the hist
  
    # correlations surrogates (null model)
    corr_surr <- mapply(cor, as.data.frame(surrogates[[(comb_two[1,i])]]), as.data.frame(surrogates[[(comb_two[2,i])]]))
    # corr_surr is a vector with 500 correlations = column by column correlation of the two matrices
  
    # binding all variables pairs in one big data.frame
    corr_pairs <- rbind.data.frame(corr_pairs, cbind.data.frame(
                       var1 = names(surrogates)[(comb_two[1,i])],
                       var2 = names(surrogates)[(comb_two[2,i])],
                       season_var1 = round(season_var1,5),
                       season_var2 = round(season_var2,5),
                       corr_data = round(corr_data,5), 
                       corr_splines = round(corr_splines,5),
                       corr_resid = round(corr_resid,5), 
                       round(t(corr_surr),5))) 
  }

  corr_pairs[1:50,1:10]
  str(corr_pairs)

  corr_pairs[,1:2] <- data.frame(lapply(corr_pairs[,1:2], as.character), stringsAsFactors=FALSE)


  if (!file.exists("output/gom_surrogates_correlation")){
    dir.create("output/gom_surrogates_correlation/") # creating a folder for CSV files
  }

  write.csv(corr_pairs, "output/gom_surrogates_correlation/correlations.csv", row.names = F)

  
  
  
  
  
# PLOTS pairs correlations ***********************************************************************************************************************

for(i in unique(c(corr_pairs$var1,corr_pairs$var2))){ # i = focus variable
  corr_subset <- corr_pairs[c(which(corr_pairs$var1 == i), which(corr_pairs$var2 == i)),]
  
  if (!file.exists(paste("output/gom_surrogates_correlation/",i,sep=""))){
    dir.create(paste("output/gom_surrogates_correlation/",i,sep="")) # creating a folder for CSV files
  }
  
  for(j in 1:nrow(corr_subset)){
    comp_var <- corr_subset[j,which(corr_subset[j,1:2]!=i)] # comparison variable
    surrgs <- data.frame(t(corr_subset[j,8:ncol(corr_subset)]), row.names = NULL)

    
    hist_corr <- ggplot(data=surrgs, aes(x = surrgs[,1]))  + 
      geom_histogram(aes(y=..density..),binwidth = 0.05, colour="grey20", fill="grey90") + # Histogram with density instead of count on y-axis
      # geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
      geom_vline(xintercept=corr_subset[j,"corr_data"],    color = "red",   size =1) +
      geom_vline(xintercept=corr_subset[j,"corr_resid"],   color = "blue",  linetype="dashed", ) +
      geom_vline(xintercept=corr_subset[j,"corr_splines"], color = "black", linetype="dashed") +
      # geom_vline(xintercept=0, color = "black") +
      labs(x = "Correlation", y = "Density") +
      theme_bw() + xlim(-1, 1) +
      theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
            strip.text = element_text(size = 20, face="bold"))
    
    pdf(file = paste("output/gom_surrogates_correlation/",i,"/",i,"-",comp_var,".pdf",sep=""), width=10, height=4, paper = "special")
      print(hist_corr)
    dev.off()  
    
  }
}


 