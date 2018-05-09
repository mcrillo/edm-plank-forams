rm(list=ls())

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./time-series-forams/")

# Libraries
# source("R/library.R")
library(R.utils) # sourceDirectory
library(lubridate) # pkg for date calculations

# Auxiliary functions
sourceDirectory("./R/aux_functions", modifiedOnly=FALSE)

# Data
gom_data <- get_gom_sst(overwrite = F)

# Differentiating data
datafd <- as.data.frame(apply(gom_data[,8:ncol(gom_data)], 2, diff))
datafd <- cbind(open=ymd(gom_data[-1,"open"]), datafd)
# Excluding first difference between samples that had a gap inbetween bigger than 10 days
days_closed = 10 
if (any(gom_data$next_open > days_closed)){
  rows <- which(gom_data$next_open > days_closed)
  datafd <- datafd[-rows,] 
}

### ORIGINAL DATA OR DIFFERENTIATED?
gom_columns <- gom_data[,-c(2:7)] # original data
# gom_columns <- datafd # first diffERENCE
gom_columns$open <- ymd(gom_columns$open) # transforming to date format

library(GGally)
ggpairs(gom_columns[,2:ncol(gom_columns)])



# STATIONARY ***********************************************************************************************************************

# Testing if time-series are stationary
library(tseries)
# ADF
adf_series <- data.frame(splines = rep(NA,ncol(splines)), raw_data = rep(NA,ncol(splines)))
for (i in 2:ncol(gom_columns)){
  adf_series[i,"splines"] <- adf.test(splines[,i])$p.value
  adf_series[i,"data"] <- adf.test(gom_columns[,i])$p.value
}
any(adf_series > 0.05)
adf_series <- cbind(names(splines),adf_series)
# KPSS
kpss_series <- data.frame(splines = rep(NA,ncol(splines)), raw_data = rep(NA,ncol(splines)))
for (i in 2:ncol(gom_columns)){
  kpss_series[i,"splines"] <- kpss.test(splines[,i])$p.value
  kpss_series[i,"data"] <- kpss.test(gom_columns[,i])$p.value
}
any(kpss_series < 0.05)
kpss_series <- cbind(names(splines),kpss_series)


# CROSS-CORRELATION (LAGS) ***********************************************************************************************************************

# Cross-Correlation
# http://www.michaeljgrogan.com/cross-correlation-r/
ccf_all <- data.frame()
comb_two <- combn(2:ncol(gom_columns),2)

for(i in 1:length(comb_two[1,])){
  cross_corr <- ccf(gom_columns[,comb_two[1,i]], gom_columns[,comb_two[2,i]], plot=FALSE)
  cross_df <- data.frame(lag = cross_corr$lag, cr_corr = cross_corr$acf, cr_corr_abs = abs(cross_corr$acf))
  
  # Regression (think about negative lag values and who is y
  lag1 <- cross_df[which.max(cross_df$cr_corr_abs),"lag"]
  lag2 <- cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),"lag"]
  
  ccf_all <- rbind(ccf_all, cbind(var1=names(gom_columns)[comb_two[1,i]],
                                  var2=names(gom_columns)[comb_two[2,i]],
                                  cross_df[which.max(cross_df$cr_corr_abs),],
                                  cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),])) 
} # for

ccf_all[,grep(names(ccf_all), pattern=c("corr"))] <- round(ccf_all[,grep(names(ccf_all), pattern=c("corr"))],2)
write.csv(ccf_all, "output/gom_ccf.csv", row.names = F)

# Lag 2 (-39)
# yplus39=tail(diffnumberpeople,26025)
# reg2<-lm(yplus39~difftemperature[1:26025])
# summary(reg2)

# reg2 <- lm(y ~ x)
# summary(reg2)$adj.r.squared
# summary(reg2)$r.squared
# summary(reg2)$coefficients[2,1] # Estimate slope regression 
# summary(reg2)$coefficients[2,4] # p-value
# library(lmtest)
# dwtest(reg)



# SPLINES & SURROGATES ***********************************************************************************************************************

### Splines
splines <- as.data.frame(SeasonalSplines_marina(DataSeries = gom_columns, 
                                                DateFormat='%d/%m/%y', # ORIGINAL: '%d/%m/%Y',
                                                SavePlots=T), stringsAsFactors = FALSE)

splines[,2:ncol(splines)]<- sapply(splines[,2:ncol(splines)], as.numeric)
splines$open <- dmy(splines$open) # transforming to date format
# data.frame(names(gom_columns), names(splines))

# Correlation of raw data and spline for each species
seasonality_species <- data.frame(species = rep(NA,ncol(splines)), corr_spline = rep(NA,ncol(splines)))
for (i in 2:ncol(splines)){
  seasonality_species[i,"species"] <- names(splines)[i]
  seasonality_species[i,"corr_spline"] <- cor(splines[,i], gom_columns[,i])
}


### Surrogate 

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


 