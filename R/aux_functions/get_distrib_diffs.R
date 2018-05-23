
# Calculates differences and sums between two consecutive time-steps of the flux data
# This is to have a distribution to sample the difference between consecutive samples (NA_resolution)

# days_closed: xcluding first differences and sums of consecutive samples that had a gap inbetween bigger than days_closed = 10 days (next_open > 10 days)


get_distrib_diffs <- function(data, trap_name, days_closed, overwrite){ # data.frame, string
  
  
  # removing NAs
  data <- data[complete.cases(data),]
       
  for(j in 1:nrow(data)){ # j = 16
    
    if(j != nrow(data)){
      data$next_open[j] <- data$open[j+1] - data$close[j]
    }else{
      data$next_open[j] <- NA
    }
  }

  # subset to just fluxes
  data_fluxes <- data[,5:(ncol(data)-1)]
  # calculating differences and sums of consecutive samples
  data_diffs <- as.data.frame(apply(data_fluxes , 2, diff))
  data_sums <- as.data.frame(apply(data_fluxes , 2, function(x) rollapply(x, 2, by = 1, sum)))
  
  # Excluding first differences and sums of consecutive samples that had a gap inbetween bigger than days_closed = 10 days (next_open > 10 days)
  if (any(data$next_open > days_closed)){
    rows <- which(data$next_open > days_closed)
    data_diffs <- data_diffs[-rows,] 
    data_sums <- data_sums[-rows,] 
  }

data_list <- list(diffs = data_diffs, sums = data_sums)
  
  if(overwrite == TRUE | !file.exists(paste("data/",trap_name,"/distrib_diffs",sep=""))){
    
    if(!file.exists(paste("data/",trap_name,"/distrib_diffs",sep=""))){ dir.create(paste("data/",trap_name,"/distrib_diffs/",sep="")) }

    for (i in 1:ncol(data_diffs)){

      savename<-paste("data/",trap_name,"/distrib_diffs/diff_",names(data_diffs)[i],".png",sep="")
      png(savename, width = 800, height = 600)
      print(ggplot(data=data_diffs, aes(data_diffs[,i])) + geom_histogram(binwidth=10) + 
        labs(x=paste(names(data_diffs)[i]," Difference between samples", sep = ""), y = "Counts"))
      dev.off()

     }
  }

 return(data_list)

}



