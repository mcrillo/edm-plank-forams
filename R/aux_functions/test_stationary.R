

test_stationary <- function(data, trap_name, overwrite){ 
  
  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/stationary_test_",trap_name,".csv",sep=""))){
    

    statio <- data.frame()
    
    # ADF
    for (i in 2:ncol(data)){
      statio[i, "variable"] <- names(data[i])
      
      # ADF
      adf_test <- adf.test(data[,i])
      statio[i,"adf_method"] <- adf_test$method
      statio[i,"adf_stats"] <- adf_test$statistic
      statio[i,"adf_lag_order"] <- adf_test$parameter
      statio[i,"alternative_hyp"] <- adf_test$alternative
      statio[i,"adf_p_value"] <- adf_test$p.value
      
      # KPSS
      kpss_test <- kpss.test(data[,i])
      statio[i,"kpss_method"] <- kpss_test$method
      statio[i,"kpss_level"] <- kpss_test$statistic
      statio[i,"kpss_trunc_lag"] <- kpss_test$parameter
      statio[i,"kpss_p_value"] <- kpss_test$p.value

    }

    statio <- statio[-1,]
    
    write.csv(statio, paste("output/",trap_name,"/stationary_test_",trap_name,".csv",sep=""), row.names = FALSE)
    
    return(statio)
    
    #paste("data/",trap_name,"/",trap_name,"_first_diff.csv",sep="")
    
    
  }else{
    
    statio <- read.csv(paste("output/",trap_name,"/stationary_test_",trap_name,".csv",sep=""), header = TRUE)
    return (statio)
  }
  
}