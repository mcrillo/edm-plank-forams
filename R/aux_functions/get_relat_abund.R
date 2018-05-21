
get_relat_abund <- function(data, trap_name, overwrite){ # data.frame, string
  
  if(overwrite == TRUE | !file.exists(paste("data/",trap_name,"/",trap_name,"_relat_abund.csv",sep=""))){
    
    # Subsetting data to just shell abundance columns
    datarel <- data[,(which(colnames(data)=="next_open")+1):ncol(data)]
    if(any(colnames(datarel)=="sst")) datarel <- datarel[,-which(colnames(datarel)=="sst")]
    
    # Calculating total shell abundance per sample
    total_abund <- data.frame(total_abund = rowSums(datarel))
    
    # Calculating relative abundances per sample
    datarel <- datarel/total_abund$total_abund
    # rowSums(datarel) # double-check
    
    datarel <- cbind(open=ymd(data[,"open"]), cum_days = data[,"cum_days"], datarel, total_abund)
    plot(datarel$total_abund ~ datarel$cum_days, type = "l")
    
    write.csv(datarel, paste("data/",trap_name,"/",trap_name,"_relat_abund.csv",sep=""), row.names = FALSE)
    return (datarel)
    
  }else{
    
    datarel <- read.csv(paste("data/",trap_name,"/",trap_name,"_relat_abund.csv",sep=""), header = TRUE)
    return (datarel)
  }
  
}
