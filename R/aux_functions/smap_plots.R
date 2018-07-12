

# data_ts: data.frame as time-series class
# emax: vector with maximum embedding dimensions to use 
# trap_name: name of sediment trap
# overwrite: if TRUE, re-runs the function, if FALSE only runs the function if plots do not exist


smap_plots <- function(data_ts, emax, trap_name, overwrite){    
  
  if(overwrite == TRUE | !file.exists(paste("output/",trap_name,"/smap_plots",sep=""))){
    
    if (!file.exists(paste("output/",trap_name,"/smap_plots",sep=""))){ dir.create(paste("output/",trap_name,"/smap_plots/",sep="")) }
    
    theta <- c()
    
    for (i in 2:ncol(data_ts)){
      X <- as.data.frame(data_ts)[,i]
      
      smap_output <- s_map(X, lib = c(1, NROW(X)), pred = c(1, NROW(X)), E = emax[i-1], silent = T)
      theta[i-1] <- which.max(smap_output$rho)
      
      savename<-paste("output/",trap_name,"/smap_plots/theta_",colnames(data_ts)[i],".png",sep="")
      png(savename, width = 800, height = 600)
      par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
      plot(x=smap_output$theta, y=smap_output$rho, type = "b", xlab = expression(paste("Theta (", sep = "", theta, ")")), 
           ylab = "Forecast Skill (rho)")
      points(y=smap_output$rho[(theta[i-1])], x=smap_output$theta[(theta[i-1])], col = "red", lwd = 4, pch = 19)
      dev.off()
    }

    print(data.frame(species = colnames(data_ts)[-1], theta))
  }
  
  
} # function
