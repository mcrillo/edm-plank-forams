# code from http://www.michaeljgrogan.com/cross-correlation-r/

cross_correlation <- function(data, trap_name, overwrite){ 
  
  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""))){
    
    
    ccf_all <- data.frame()
    
    # combination of all pairs of variables to test for cross-correlation and lags
    comb_two <- combn(2:ncol(data),2)
    
    for(i in 1:length(comb_two[1,])){ # for each pair of variables
      
      cross_corr <- ccf(data[,comb_two[1,i]], data[,comb_two[2,i]], plot=FALSE)
      cross_df <- data.frame(lag = cross_corr$lag, cr_corr = cross_corr$acf, cr_corr_abs = abs(cross_corr$acf))
      
      ccf_all <- rbind(ccf_all, cbind(var1=names(data)[comb_two[1,i]],
                                      var2=names(data)[comb_two[2,i]],
                                      cross_df[which.max(cross_df$cr_corr_abs),],
                                      cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),],
                                      cross_df[which(cross_df$lag == 0),]
                                      ) ) 
      
      
    } # for
    
      # Regression (think about negative lag values and who is y
      
      #lag1 <- cross_df[which.max(cross_df$cr_corr_abs),"lag"]
      #lag2 <- cross_df[which.max(cross_df$cr_corr_abs[cross_df$cr_corr_abs!=max(cross_df$cr_corr_abs)]),"lag"]
      
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

    ccf_all[,grep(names(ccf_all), pattern=c("corr"))] <- round(ccf_all[,grep(names(ccf_all), pattern=c("corr"))],2)

    write.csv(ccf_all, paste("output/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""), row.names = FALSE)
    return(ccf_all)
    
  }else{
    
    ccf_all <- read.csv(paste("output/",trap_name,"/cross_correlation_",trap_name,".csv",sep=""), header = TRUE)
    return (ccf_all)
  }
  
}