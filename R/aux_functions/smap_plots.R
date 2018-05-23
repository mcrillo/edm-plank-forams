

# data_ts: data.frame as time-series class
# emax: vector with maximum embedding dimensions to use 
# trap_name: name of sediment trap
# overwrite: if TRUE, re-runs the function, if FALSE only runs the function if plots do not exist


smap_plots <- function(data_ts, emax, trap_name, overwrite){    
  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/smap_plots",sep=""))){
    
    if (!file.exists(paste("output/",trap_name,"/smap_plots",sep=""))){ dir.create(paste("output/",trap_name,"/smap_plots/",sep="")) }
    
    for (i in 2:ncol(data_ts)){
      X <- as.data.frame(data_ts)[,i]
      
      smap_output <- s_map(X, lib = c(1, NROW(X)), pred = c(1, NROW(X)), E = emax[i-1], silent = T)
      
      savename<-paste("output/",trap_name,"/smap_plots/theta_",colnames(data_ts)[i],".png",sep="")
      png(savename, width = 800, height = 600)
      par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
      plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", 
           ylab = "Forecast Skill (rho)")
      dev.off()
    }

  }
  
  
} # function
