#1.Libraries====================================================================
library(ncdf4)
library(ff)
library(ffbase)
library(reshape2)
library(tidyverse)
library(magrittr)
library(stringr)		
library(REdaS)

years<- seq(2006,2049,1)

source("/lustre/scratch2/ws/1/ahho623a-FRM_PDFs_project/CORDEX/RCP_85/proj_CORDEX_RCP85.R")

for(i in 2:length(years)){
	
setwd(paste0("/lustre/scratch2/ws/1/ahho623a-FRM_PDFs_project/CORDEX/RCP_85/",years[i]))
	
getwd()

nc.files<-list.files(pattern = ".nc$")

nc.files<- nc.files[-1]

nc.files

system("mkdir projected")


for (i in 1:length(nc.files)) {
	
	ta850 <-grep(pattern = 'ta850', x = nc.files,value = TRUE)
	
	print(i)
	
	proj_read_toff(netcdf_file=nc.files[i],
		                refrence_file= ta850)
	
	#rm(ncin)				          # remove dummy variables 
	#rm(med)
	gc()
	
}
	
}