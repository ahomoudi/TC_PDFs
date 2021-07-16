
library(plyr)
library(ff)
library(ffbase)
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(lubridate)

library(sp)
library(maptools)
library(dplyr)
library(sf)
library(leaflet)

library(mapview)

# Read data ---------------------------------------------------------------


ERA5<- read.csv("/media/ahmed/Volume/TC_FRM/IBTrACS.NI.list.v04r00.lines/points_processed.csv")



all_content = readLines("/media/ahmed/Volume/TC_FRM/R/ibtracs.NI.list.v04r00.csv")

skip_second = all_content[-2]

IBTrACS <- read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)

rm(all_content,skip_second)

# filter data 
IBTrACS<-IBTrACS[which(IBTrACS$SEASON>=1990 & IBTrACS$SEASON<=2019),]

IBTrACS<-IBTrACS[which(IBTrACS$SUBBASIN =="AS"),]

# Read function -----------------------------------------------------------

source("/media/ahmed/Daten/GCMs/ITCZ_analysis/extra_function.R")
source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/points_to_line.R")


Points_df<- ERA5%>% select("time" ,"Season","Storm","Track", "MTS")

count_df<-Points_df[Points_df$Track==1,]%>% 
	group_by(Season)%>%
	count(MTS)

ggplot(count_df, aes(Season, y = n,fill= MTS)) +
	geom_col()+
	scale_x_discrete() 

