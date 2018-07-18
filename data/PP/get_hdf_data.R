
rm(list=ls())

#setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
#setwd("./time-series-forams/data/PP/") 
setwd("/Users/marinacostarillo/Desktop")

# Libraries
library(R.utils)
library(geosphere)

# Auxiliary functions
find_neighbours <- function(point, findin, distance) { # vector, data.frame, numeric
  ## point: vector of two numbers (longitude, latitude) 
  ## findin: a matrix of 2 columns (first one is longitude, second is latitude) 
  ## distance: if 0 finds nearest neighbour, if a positive value (in meters) finds neighbours within a circle with the value radius
  
  ## Returns a data.frame: 
  ## "row" = the row number of the neighbour in data.frame
  ## "distance" = the distance between the points
  
  dist_data <- apply(findin, 1, function(x) distCosine(point, x)) #Matrix
  
  if(distance>0) { # find neighbours within radius of distance
    neighb <- data.frame(row_findin = which(dist_data<=distance), distance = dist_data[which(dist_data<=distance)])
    if(length(neighb[,1])==0) distance = 0 
  }  
  
  if(distance==0) { # find nearest neighbour
    neighb <- data.frame(row_findin = which.min(dist_data), distance = min(dist_data))
  }
  
  return(neighb)
  
}

files <- list.files(pattern = ".gz")
i = 1
gunzip(files[i])

# vgpm.2008097.all.xyz
data_xyz <- read.delim(file = "vgpm.2008097.all.xyz", header = TRUE, sep = " ", dec = ".")
data_xyz <- data_xyz[,1:3]

rows_xyz <- find_neighbours(point = c(-90.3,27.5), findin = data_xyz[,c("lon","lat")], distance = 0)

require(raster) 
d <- raster(somefile) 

d <- rasterFromXYZ(somefile) 

require(rgdal) 
d <- readGDAL(somefile) 

# cbpm.2010097.hdf
gdalinfo("cbpm.2010097.hdf")


### rhdf5
https://stackoverflow.com/questions/15974643/how-to-deal-with-hdf5-files-in-r
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
library(rhdf5)
h5ls("cbpm.2010097.hdf") # does not work! is this hdf file what? 4, 5 ?


### hdf5r
# install.packages("hdf5r") 
# Problem:
# hdf5r is only available in source form, and may need compilation of C/C++/Fortran: 'hdf5r' 
# compilation failed for package ‘hdf5r’
# brew no terminal: Warning: hdf5 dependency gcc was built with a different C++ standard library (libstdc++ from clang). This may cause problems at runtime.
# Apple is moving away from libstdc++, the system-provided libstdc++ has not been updated since OS X 10.7 (Lion) and it lacks any C++11 facility.
# http://ease-the-computation.haunted.io/compile-clang-against-libstdc-with-c11-support-on-a-mac/

### xyz files: http://orca.science.oregonstate.edu/2160.by.4320.8day.xyz.vgpm.m.chl.m.sst.php
# week 2008097 http://orca.science.oregonstate.edu/data/2x4/8day/vgpm.r2018.m.chl.m.sst/xyz/vgpm.2008097.all.xyz.gz
