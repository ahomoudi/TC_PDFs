
library(plyr)
library(ff)
library(ffbase)
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(lubridate)
#1.Reading files########################################################################

#'* ATTENTION CHANGE  WORKING DIRCTORY *

setwd(paste0("/media/ahmed/Volume/TC_FRM/Results_ERA5"))

# Create list of text files
csv_files <- list.files(pattern="*.csv")

# Read the files in, assuming comma separator
csv_files_df <- lapply(csv_files, function(x) {
	y<-unlist(str_split(x, "[.]"))[1]
	assign(y,
	       as.ffdf(fread(x, stringsAsFactors = TRUE)))})

# Combine them
combined_df<-plyr::ldply(csv_files_df,data.frame)

# convert pr from "m" to "mm"
combined_df$prmax_grd<-combined_df$prmax_grd*1000

combined_df$prsum_grd<-combined_df$prsum_grd*1000

hist(combined_df$prsum_grd)

#remove points where no rainfall is recorded
combined_df<-combined_df[-which(combined_df$prsum_grd<4.9),]

hist(combined_df$prsum_grd)

# add season to df 
combined_df$Season<-year(combined_df$time)

# write filtered points
write.csv(combined_df,paste0("/media/ahmed/Volume/TC_FRM/Results_ERA5/Results/All.csv"))

combined_df$time<-as.POSIXct(strptime(as.character(combined_df$time), "%Y-%m-%d %H:%M:%S"))


filter_df<-data.frame(No= seq(1,nrow(combined_df),1))

filter_df$Season<-combined_df$Season

# filter according to time 
for( i in 1:nrow(combined_df)){
	if( i ==1 ){
		
	xt<-difftime(combined_df$time[i+1],combined_df$time[i],units = "days")
	
	xd<-sqrt((combined_df$LAT[i+1]-combined_df$LAT[i])**2 +
	        	(combined_df$LON[i+1]-combined_df$LON[i])**2)
	
	if(as.numeric(xd)<5 &as.numeric(xt)<1){filter_df$DT[i]<-2}else{filter_df$DT[i]<-0}
	
	
	}
	if( i == nrow(combined_df) ){
	xt<-difftime(combined_df$time[i],combined_df$time[i-1],units = "days")
	
	xd<-sqrt((combined_df$LAT[i]-combined_df$LAT[i-1])**2 +
	        	(combined_df$LON[i]-combined_df$LON[i-1])**2)
	
	if(as.numeric(xt)<1& as.numeric(xd)<5 ){filter_df$DT[i]<-1}else{filter_df$DT[i]<-0}
		
	}
	if(i !=1 & i!=nrow(combined_df)){
		
	xt<-difftime(combined_df$time[i],combined_df$time[i-1],units = "days")
	yt<-difftime(combined_df$time[i+1],combined_df$time[i],units = "days")
	
	xd<-sqrt((combined_df$LAT[i]-combined_df$LAT[i-1])**2 +(combined_df$LON[i]-combined_df$LON[i-1])**2)
	yd<-sqrt((combined_df$LAT[i+1]-combined_df$LAT[i])**2 +(combined_df$LON[i+1]-combined_df$LON[i])**2)
	
	if(as.numeric(xt)<1& as.numeric(yt)<1& as.numeric(xd)<5 & as.numeric(yd)<5){filter_df$DT[i]<-3}
	else if(as.numeric(xt)>1& as.numeric(yt)<1& as.numeric(xd)>5 & as.numeric(yd)<5){filter_df$DT[i]<-2}
	else if(as.numeric(xt)<1& as.numeric(yt)>1& as.numeric(xd)<5 & as.numeric(yd)>5){filter_df$DT[i]<-1}
	else {filter_df$DT[i]<-0}
	}
}


filtered_points_df<-combined_df[-which(filter_df$DT==0),]

# write extra filtered points

write.csv(filtered_points_df,paste0("/media/ahmed/Volume/TC_FRM/Results_ERA5/Results/All_Filtered.csv"))

filter_df2<-filter_df[-which(filter_df$DT==0 ),]

hist(filter_df2$DT)

Start_point<-which(filter_df2$DT==2)

blabla<-length(Start_point)-1

filter_df2$Storm<-NA;filter_df2$Track<-NA

for ( i in 1:blabla){
	
storm<-filter_df2[Start_point[i]:(Start_point[i+1]-1),]

x<-nrow(storm)

y<-x-1

#if(storm$TG[1]==2 & storm$TG[x]==1 & all(storm$TG[2:y]==3) & all(storm$DG>=2)){
if(storm$DT[1]==2 & all(storm$DT[2:y]==3)){
	
	filter_df2$Storm[Start_point[i]:(Start_point[i+1]-1)]<-paste0("Storm_",i)
	
	filter_df2$Track[Start_point[i]:(Start_point[i+1]-1)]<-seq(1,x,1)
	
}}

# add track information to data frame 
filtered_points_df$Storm<-filter_df2$Storm

filtered_points_df$Track<-filter_df2$Track

library(sp)
library(maptools)
library(dplyr)

source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/points_to_line.R")

Points_df<- filtered_points_df%>% select("LAT","LON", "Storm","Track")

colnames(Points_df)<-c("lat","long", "Storm","Track")

write.csv(filtered_points_df,paste0("/media/ahmed/Volume/TC_FRM/Results_ERA5/Results/Points_df.csv"))

unique(Points_df$Storm)


library(leaflet)

v_lines <- points_to_line(data = Points_df, 
	      long = "long", 
	      lat = "lat", 
	      id_field = "Storm", 
	      sort_field = "Track")

m<-leaflet(data = v_lines) %>%
	addTiles() %>%
	addPolylines()

library(mapview)
mapshot(m, url = paste0("/media/ahmed/Volume/TC_FRM/Results_ERA5/Results/tracks.html"))
















library(sf)
st_write(v_lines%>% st_as_sf(),                                                  #write polygons to shape file
         paste0(getwd(),"Results/v_lines.shp"),
         overwrite=T) 
