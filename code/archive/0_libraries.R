###  ### ### ### ### ### ### ### ### ### ### ### ###             
### BACI Analysis of MPA effects on kelp = Packages
###### ### ### ### ### ### ### ### ### #### ### ###    

library(tidyverse);
library(car)
library(cowplot);
library(here);
library(wesanderson);
library(tidyr)
library(readr);
library(MuMIn);
library(here);
library(lubridate);
library(data.table);
library(ggpubr);
library(ncdf4); # package for netcdf manipulation
library(raster); # package for raster manipulation
library(rgdal); # package for geospatial analysis
library(ggplot2); # package for plotting
library(pacman);
library(sf);
library(data.table);
library(mapview);
library(fields);
library(mgcv);
library(olsrr);
library(emmeans);
library(stats);
library(Rmisc);
library(gtsummary);
library(broom);
library(gt);
library(lme4);
library(lmtest);
library(lmerTest);
library(dplyr);
library(car);
library(vegan);
require(minpack.lm); # Fitting non-linear models
require(nls2); # Fitting non-linear models
require(AICcmodavg); # calculate second order AIC (AICc)
library(remef); # to extract leverage pairs
library(nlstools);# to run autocorrelation test on resids of nls2 model
library(investr); # to extrace predictions from nls model
library(rempsyc); # exporting tables
library(flextable)
library(yarrr) #for color transparency
library(formattable);
library(googlesheets4) #For importing PISCO data
library(plyr)
library(metafor)
pal <- wes_palette("Zissou1", 25, type = "continuous")
