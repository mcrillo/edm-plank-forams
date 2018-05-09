

surrogates_splines <- function(data, splines, trap_name, nreps, overwrite){     

  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/correlation_surrogates_",trap_name,".csv",sep=""))){
    
    # (1) # Randomization of residuals of (data - splines)
    # data.frame(data = data[,i], spline = splines[,i], resid = data[,i] - splines[,i])
    random_resid <- list()
    random_resid <- lapply(2:ncol(splines), function(i) replicate(nreps, sample(c(data[,i] - splines[,i]), size = length(splines[,i]), replace = FALSE, prob = NULL)))
    names(random_resid) <- names(splines)[2:ncol(splines)]
    # str(random_resid) # Each element of the list is a matrix with 500 randomized residual series
    # var(random_resid[[i]][,90]) # double check if variance is always the same, i.e., if vectors are all identical just with values in diferent order (randomized)
    
    
    # (2) # Generating surrogates for each species: adding random residuals to splines (i.e., seasonal time-series)
    # data.frame(resid = random_resid[[i]][,500], spline = splines[,i],  surrogate = random_resid[[i]][,500]+ splines[,i])
    surrogates <- lapply(1:length(random_resid), function(i){random_resid[[i]] + splines[,i+1]}) # i+1 because first column of splines is "open"
    names(surrogates) <- names(random_resid)
    # str(surrogates) # Each element of the list is a matrix with 500 null series (= spline + random residual)
    # all(round(surrogates[[i]][,500] - random_resid[[i]][,500] - splines[i+1],4)==0) # double check: this must be all 0 for any value of i in [1,15] and j[1,500]
    
    # Calculate correlations column by columns of each variable pair
    corr_pairs <- data.frame()
    comb_two <- combn(1:length(surrogates),2)
    
    for(i in 1:length(comb_two[1,])){
      
      col_var1 <- which(names(data)==names(surrogates)[(comb_two[1,i])])
      col_var2 <- which(names(data)==names(surrogates)[(comb_two[2,i])])
      
      # correlation of each variable with its own spline (measure of seasonality of the time-series)
      season_var1 <- cor(data[,col_var1], splines[,col_var1])
      season_var2 <-  cor(data[,col_var2], splines[,col_var2])
      
      # correlations
      corr_series <- cor(data[,col_var1], data[,col_var2])
      corr_splines <- cor(splines[,col_var1], splines[,col_var2]) 
      corr_resid <- cor( (data[,col_var1]-splines[,col_var1]), (data[,col_var2]-splines[,col_var2]) ) 
      
      # correlations surrogates (null model)
      corr_surr <- mapply(cor, as.data.frame(surrogates[[(comb_two[1,i])]]), as.data.frame(surrogates[[(comb_two[2,i])]]))
      # correlates each column of the matrices for each species (500 columns each, null series)
      # corr_surr is a vector with 500 correlations = column by column correlation of the two matrices
      
      # significance
      probs <- quantile(corr_surr, probs = c(0.025, 0.975)) 
      if(corr_series < probs[1] | corr_series > probs[2]) {
        p_series <- c("signif")
      }else{
        p_series <- c("non")
      }
      
      # binding all variables pairs in one big data.frame
      corr_pairs <- rbind.data.frame(corr_pairs, cbind.data.frame(
        var1 = names(surrogates)[(comb_two[1,i])],
        var2 = names(surrogates)[(comb_two[2,i])],
        season_var1 = round(season_var1,5),
        season_var2 = round(season_var2,5),
        corr_series = round(corr_series,5), 
        p_series = p_series , 
        corr_splines = round(corr_splines,5),
        # p_splines = NA,
        corr_resid = round(corr_resid,5), 
        # p_resid = NA,
        round(t(corr_surr),5))) 
    }
    
    # corr_pairs[1:50,1:10]
    # str(corr_pairs)
    
    corr_pairs[,1:2] <- data.frame(lapply(corr_pairs[,1:2], as.character), stringsAsFactors=FALSE)
    
    write.csv(corr_pairs, paste("output/",trap_name,"/correlation_surrogates_",trap_name,".csv",sep=""), row.names = F)
    return(corr_pairs)
    
  }else{
    
    corr_pairs <- read.csv(paste("output/",trap_name,"/correlation_surrogates_",trap_name,".csv",sep=""), header = TRUE)
    return (corr_pairs)
  }
  
}