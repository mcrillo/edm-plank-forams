
# days_closed: xcluding first differences and sums of consecutive samples that had a gap inbetween bigger than days_closed = 10 days (next_open > 10 days)


get_first_diff <- function(data, trap_name, days_closed, overwrite){ # data.frame, string
  
  if(overwrite == TRUE | !file.exists(paste("data/",trap_name,"/",trap_name,"_first_diff.csv",sep=""))){
      
      # Differentiating data
      datafd <- as.data.frame(apply(data[,(which(colnames(data)=="next_open")+1):ncol(data)], 2, diff))
      datafd <- cbind(open=ymd(data[-1,"open"]), datafd)
      
      # Excluding first difference between samples that had a gap inbetween bigger than days_closed = 10 days (next_open > 10 days)
      if (any(data$next_open > days_closed)){
        rows <- which(data$next_open > days_closed)
        datafd <- datafd[-rows,] 
      }
      
      write.csv(datafd, paste("data/",trap_name,"/",trap_name,"_first_diff.csv",sep=""), row.names = FALSE)
      return (datafd)
  
   }else{
  
       datafd <- read.csv(paste("data/",trap_name,"/",trap_name,"_first_diff.csv",sep=""), header = TRUE)
       return (datafd)
   }
  
}

