
# data_ts: data.frame as time-series class
# emax: vector with maximum embedding dimensions to use 
# trap_name: name of sediment trap
# overwrite: if TRUE, re-runs the function, if FALSE only runs the function if plots do not exist


simplex_plot <- function(data_ts, emax, trap_name, overwrite){    
  
  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/simplex_pred_plots",sep=""))){
    
    if (!file.exists(paste("output/",trap_name,"/simplex_pred_plots",sep=""))){ dir.create(paste("output/",trap_name,"/simplex_pred_plots/",sep="")) }
    
  
 rho_df <- data.frame()

 for (i in 2:ncol(data_ts)){
  
  X <- as.data.frame(data_ts)[,i]
  pred <- simplex(X, lib = c(1, NROW(X)), pred = c(1, NROW(X)), E = emax[i-1], tp = 1:20, stats_only = FALSE, silent = T)
  fits <- pred$model_output[[1]]
  # head(fits)
  
  rho_df[i,"variable"] <- colnames(data_ts)[i]
  rho_df[i,"rho"] <- round(pred$rho[1],3)
  rho_df[i,"emax"] <- emax[i-1]
  
  savename1<-paste("output/",trap_name,"/simplex_pred_plots/",colnames(data_ts)[i], "_pred.png",sep="")
  png(savename1, width = 800, height = 500)
  plot(pred ~ time, data = fits, type = "l", col = "blue", lwd=3,
       xlab="Time", ylab=colnames(data_ts)[i], ylim=range(fits[,2:3],na.rm = T))
  lines(obs ~ time, data = fits, col=grey.colors(1, alpha=0.25), lwd = 6)
  legend("topright", c("Observed", "Predicted"), lty=1, lwd=c(6,3),
         col=c(grey.colors(1, alpha=0.25), "blue"),bty="n")
  dev.off()
  
  savename2<-paste("output/",trap_name,"/simplex_pred_plots/decay_",colnames(data_ts)[i], ".png",sep="")
  png(savename2, width = 800, height = 600)
  plot(rho ~ tp, data=pred,
       type = "b",
       xlab = "Time to prediction",
       ylab = expression(paste("Forecast skill (",rho,")",sep="")))
  dev.off()
  
}
rho_df <- rho_df[-1,]

write.csv(rho_df, paste("output/",trap_name,"/simplex_rho_",trap_name,".csv",sep=""), row.names = F)
return(rho_df)


  }else{
    
    rho_df <- read.csv(paste("output/",trap_name,"/simplex_rho_",trap_name,".csv",sep=""), header = TRUE, stringsAsFactors = FALSE)
    return (rho_df)
  }
  
  
} # function

