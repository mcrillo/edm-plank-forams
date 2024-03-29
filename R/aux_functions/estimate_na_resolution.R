

estimate_na_resolution <- function(data, list, overwrite){

if(overwrite == TRUE | !file.exists("data/GOM/GOM_edm_use.csv")){
    
  # time-steps that need to be estimated
  data_na_resolution <- data[duplicated(data$open),]

  # samples averaged throughout the 14 days
  samples_resolution <- as.numeric(as.character(rownames(data_na_resolution)))-1
  data_samples <- data[samples_resolution,]

  for (i in 5:(ncol(data_na_resolution)-1)){ # for each species / column 
    # print(names(data_na_resolution)[i])
    for (j in 1 :nrow(data_na_resolution)){
      
      if(data_samples[j,i]!=0){
      # if NA_resolution did not existed, and the two weeks were samples normally, here is what you would have expected as the sum of shell of the two consecutives samples
        total_shells <- 2*data_samples[j,i]
        # print(total_shells)
        rows_to_sample <- which(list$sums[,which(names(list$sums)==names(data_na_resolution)[i])]<=total_shells & list$sums[,which(names(list$sums)==names(data_na_resolution)[i])]!=0)

        diff_distrib <- list$diffs[rows_to_sample,which(names(list$sums) == names(data_na_resolution)[i])]
      
        # repeat until both samplings are positive values (= shell flux always positive)
        if(any(diff_distrib<data_samples[j,i])){
          set.seed(42) #for reproducibility
          while(TRUE){
            peteleco <- sample(diff_distrib,1) 
            sample1 <- data_samples[j,i]+peteleco
            sample2 <- data_samples[j,i]-peteleco
            if(sample1>=0 && sample2>=0 && sample1!=sample2) break()
          }
        }else{
          peteleco <- runif(1, min=0, max=data_samples[j,i])
          sample1 <- data_samples[j,i]+peteleco
          sample2 <- data_samples[j,i]-peteleco
        }

        time_step1 <- data_samples[j,"time_step"]
        time_step2 <- data_na_resolution[j,"time_step"]
      
        if(sample1+sample2!=total_shells){"Check code 1"} # MEDIA!
    
        data[which(data$time_step==time_step1),i] <- sample1
        data[which(data$time_step==time_step2),i] <- sample2
    
        rm(rows_to_sample,diff_distrib,peteleco,sample1,sample2)
    
      }else{ # if 0 shells were sampled
        time_step2 <- data_na_resolution[j,"time_step"]
        data[which(data$time_step==time_step2),i] <- data_samples[j,i]
      }
    }
  
    if(any(data[,i]<0, na.rm = T)){"Check code 2"}
  
  } # species
  
  write.csv(data, "data/GOM/GOM_edm_use.csv", row.names = FALSE)
  return (data)
  
}else{
    
    data <- read.csv("data/GOM/GOM_edm_use.csv", header = TRUE)
    return (data)
}

} # function




