# get Daily_SST_GOM

  rm(list=ls())

  library(geosphere)
  library(lubridate)
  library(ncdf4)
  # Reynolds, R. W., T. M. Smith, C. Liu, D. B. Chelton, K. S. Casey and M. G. Schlax, 2007:Daily High-resolution Blended Analyses for sea surface temperature. J. Climate, 20, 5473-5496. | Reynolds, R. W., 2009:What's New in Version 2. http://www.ncdc.noaa.gov/sst/papers/oisst_daily_v02r00_version2-features.pdf
  setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
  setwd("./time-series-forams/data/OISST/")  
  oisst_files <- list.files("data_oisst/")
  
  ### Finding nearest point to GoM sediment trap
  # GoM latN: 27.5	longE: -90.3 = 269.7
  points_oisst <- data.frame(lon=c(269.625,269.625,269.875,269.875), 
                             lat=c(27.375,27.625,27.375,27.625))
  min(c(
    distCosine(p1=c(269.7,27.5),p2=points_oisst[1,]), # 15764.86
    distCosine(p1=c(269.7,27.5),p2=points_oisst[2,]), # 15760.91 nearest
    distCosine(p1=c(269.7,27.5),p2=points_oisst[3,]), # 22193.58
    distCosine(p1=c(269.7,27.5),p2=points_oisst[4,]) # 22178.29
  ))
  
  which(data$lat == points_oisst[p,"lat"])
  which(data$lon == points_oisst[p,"lon"])
  
  # start date of oisst (metadata)
  start <- ymd(19780101)

  gom_daily_sst <- data.frame(temp_K = double(), day = double())

  for (i in oisst_files){ # i = oisst_files[1]
    
      print(i)
      
      # open a NetCDF file
      nc <- nc_open(paste("data/",i, sep=""))
      # nc
    
      data <- list()
      data$lat <- ncvar_get(nc, "lat") # units: degrees_north
      data$lon <- ncvar_get(nc, "lon") # units: degrees_east
      data$time <- ncvar_get(nc, "time") # units: days since 1978-01-01 00:00:00
      data$tos <- ncvar_get(nc, "tos") # units: K. Comment: converted to units in Kelvin by subtracting 273.15 from original values in Celsius.
    
      # str(data) # tos[lon,lat,time]   (Chunking: [1440,720,1])  
      sst <- data$tos[which(data$lon == 269.625),which(data$lat == 27.625),] # GoM (distance 15760.91 meters)
    
      # Finding real days within data$time
      time_min <- start + min(data$time)
      time_max <- start + max(data$time)
      days <- seq.Date(from=time_min, to=time_max, by="day")
    
      gom_year_sst <- data.frame(temp_K = sst, day = days)
      write.csv(gom_year_sst,paste("gom_year_sst/",i,".csv", sep=""), row.names = FALSE)
    
      gom_daily_sst <- rbind(gom_daily_sst, gom_year_sst)
      
      rm(nc)
      rm(data)
      rm(gom_year_sst)

  }
  # length(gom_daily_sst$day) # check : 365*6.5
  write.csv(gom_daily_sst,"gom_daily_sst.csv", row.names = FALSE)
  # str(gom_daily_sst)

  gom_data <- read.csv("GOM.csv", header = TRUE)
  
  for(j in 1:length(gom_data[,1])){ # j = 1
  
    open <- ymd(gom_data$open)[j]
    close <- ymd(gom_data$close)[j]
    year_days <- ymd(gom_daily_sst$day)
    
    rows <- which(ymd(gom_daily_sst$day) >= open & ymd(gom_daily_sst$day) <= close)
    gom_data[j,"temp_K"] <- mean(gom_daily_sst[rows,"temp_K"])
    
  }
  gom_data[,"temp_oC"] <- gom_data[,"temp_K"] - 273.15
  
  write.csv(gom_data,"GOM_sst.csv", row.names = FALSE)
  
    