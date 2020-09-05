#0###############################################################################
# This code is written to calculate PDFs and statistical properties of the chosen 
#variables (vorticity, surface wind speed, and vertically integrated anomalous 
# temperature) according to :
#         Camargo, S. J., & Zebiak, S. E. (2002). Improving the Detection and 
#         Tracking of Tropical Cyclones in Atmospheric General Circulation 
#             Models. Weather and   Forecasting, 17(6), 1152–1162.
#         https://doi.org/10.1175/1520-0434(2002)017<1152:ITDATO>2.0.CO;2
#                  Written by :Ahmed Homoudi
#                  September 2020
#1.Libraries====================================================================
library(ncdf4)
library(ff)
library(ffbase)
library(reshape2)
library(tidyverse)
library(magrittr)
library(Rmpi)
library(parallel)

workdir<-system("pwd",intern = TRUE)
setwd(workdir)
