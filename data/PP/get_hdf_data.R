
rm(list=ls())

library(R.utils)

setwd("/Users/marinacostarillo/Google Drive/PhD/projects")
setwd("./time-series-forams/data/PP/")  
files <- list.files(pattern = ".gz")

### hdf files 
# gdalinfo("cbpm_2008001.hdf")
install.packages("hdf5r") 
# Problem:
# hdf5r is only available in source form, and may need compilation of C/C++/Fortran: 'hdf5r' 
# compilation failed for package ‘hdf5r’
# brew no terminal: Warning: hdf5 dependency gcc was built with a different C++ standard library (libstdc++ from clang). This may cause problems at runtime.
# Apple is moving away from libstdc++, the system-provided libstdc++ has not been updated since OS X 10.7 (Lion) and it lacks any C++11 facility.
# http://ease-the-computation.haunted.io/compile-clang-against-libstdc-with-c11-support-on-a-mac/

### xyz files: http://orca.science.oregonstate.edu/2160.by.4320.8day.xyz.vgpm.m.chl.m.sst.php
# week 2008097 http://orca.science.oregonstate.edu/data/2x4/8day/vgpm.r2018.m.chl.m.sst/xyz/vgpm.2008097.all.xyz.gz
i = 1
gunzip(files[i])
vgpm.2008097.all.xyz.gz
