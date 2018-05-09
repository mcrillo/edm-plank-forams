

        
get_gom_sst <- function(overwrite){ # data.frame, string
  
  if(overwrite == TRUE | !file.exists("data/GOM/GOM_original_sst.csv")){
    
      data <- read.csv("data/GOM/GOM_original.csv", header = TRUE)
    
     # Number of species in the trap
      col_ssp <- c(3:ncol(data))
      
      
    # Getting day of the year (doy) for date in the middle of sampling interval
      data$open <- dmy(data$open)
      data$close <- dmy(data$close)
      
      data$middle <- data$open + (data$close-data$open)/2
      data$doy <- yday(data$middle)
      data$cum_days <- NA
      data$resolution <- data$close-data$open
      data$next_open <- NA
      data$sst <- NA
      
    # Getting mean daily temperature for sampling interval
      data_daily_sst <- read.csv("data/GOM/Daily_SST_GOM.csv", header = TRUE)
      data_daily_sst$temp_C <- data_daily_sst$temp_K - 273.15
      
      for(j in 1:nrow(data)){ # j = 2
        
        open <- ymd(data$open)[j]
        close <- ymd(data$close)[j]
        
        data$cum_days[j] <- data$close[j] - data$open[1]
        
        
        if(j != nrow(data)){
          data$next_open[j] <- data$open[j+1] - data$close[j]
        }else{
          data$next_open[j] <- NA
        }
        
        
        rows <- which(ymd(data_daily_sst$day) >= open & ymd(data_daily_sst$day) <= close)
        
        data[j,"sst"] <- mean(data_daily_sst[rows,"temp_C"]) 
      }
      
      # Re-order columns
      order <- c(1,2,((ncol(data)-5):ncol(data)), col_ssp)
      data <- data[,order]
      
      # Save new data
      write.csv(data, "data/GOM/GOM_original_sst.csv", row.names = FALSE)
      
      return (data)
      
  }else{
    
    data <- read.csv("data/GOM/GOM_original_sst.csv", header = TRUE)
    
    return (data)
    
    
  }
  
}