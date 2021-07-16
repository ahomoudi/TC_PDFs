library(ff)
library(ffbase)
library(data.table)
library(stringr)
library(reshape2)
library(viridis)
library(tidyverse)

# Put in your actual path where the text files are saved

workdir<-system("pwd",intern = TRUE)
setwd(workdir)

# Create list of text files
txt_files <- list.files(pattern="*.txt")

#create csv files names 
csv_files <- lapply(txt_files, function(x){
	paste0(unlist(str_split(x, "[.]"))[1],".csv")})
#transform txt to csv 
command<-paste0("cat ",txt_files," | tr -s '[:blank:]' ',' > ",csv_files)

lapply(command, function(x){
	system.time(system(x))
})
# Create list of text files
csv_files <- list.files(pattern="*.csv")

csv_files_df <- list.files(pattern = "*.csv") %>% 
	map_df(~fread(.))

combined_df<-csv_files_df

nop<- which(combined_df$warm_cores[]==0 |combined_df$min_press[]==0)

combined_df <-combined_df[-nop,]



fwrite(combined_df, file = "all_filtered_CORDEX.csv", sep = ",")

# on local machine 

workdir<-system("pwd",intern = TRUE)
setwd(workdir)

csv_file_final<-read.csv("all_filtered_CORDEX.csv",sep = ",",header = T)

# remove points not in the tropical cyclone season
library(PCICt)
library(lubridate)
time.origin <- PCICt::as.PCICt.default("1949-12-01 00:00:00", 
	                   cal = "365_day")
csv_file_final$month<- month(PCICt::as.PCICt(csv_file_final$TIME*3600*24, 
		          origin = time.origin,cal = "365_day"))
nop<-which(csv_file_final$month!=c(5,6,9,10,11))

plot(csv_file_final$month)

Final_filtered<-csv_file_final[-nop,]

plot(Final_filtered$month)

statistic_df<- data.frame(row.names = colnames(Final_filtered))

for(i in 1: ncol(Final_filtered)){
	statistic_df$mean[i]<- mean(as.numeric(unlist(Final_filtered[i][])))
	statistic_df$std[i] <- sd(as.numeric(unlist(Final_filtered[i][])))
	statistic_df$med[i] <- median(as.numeric(unlist(Final_filtered[i][])))
}


#





























# Read the files in, assuming comma separator
#csv_files_df <- lapply(csv_files, function(x) {
#y<-unlist(str_split(x, "[.]"))[1]
#	assign(y,
#   as.ffdf(fread(x)))})

# Combine them
#combined_df <- do.call("ffdfappend", lapply(csv_files_df, as.ffdf))

##combined_df <- do.call("rbind", lapply(csv_files_df, as.data.frame)) 

#statistic_df<- data.frame(row.names = colnames(combined_df))

#for(i in 1: ncol(combined_df)){
#	statistic_df$mean[i]<- mean(as.numeric(unlist(combined_df[i][])))
#	statistic_df$std[i] <- sd(as.numeric(unlist(combined_df[i][])))
#	statistic_df$med[i] <- median(as.numeric(unlist(combined_df[i][])))
#}

#combined_df$vor_std<-as.ff(combined_df$vorticity/statistic_df["vorticity","std"])

#combined_df$T_std<- as.ff(combined_df[,"anomalous"]/statistic_df["anomalous","std"])

#library(ggplot2)

#plotting_df<-fortify(combined_df[,-c(1,2,3,4,5)])

#ggplot(stack(plotting_df[,c("temp300_anomaly","temp500_anomaly",
#	        "temp700_anomaly","temp850_anomaly")]), 
#       aes(x = ind, y = values)) +geom_violin()

#ggplot(stack(plotting_df[,c("temp300","temp500",
#	        "temp700","temp850")]), 
#       aes(x = ind, y = values)) +geom_violin()

#ggplot(plotting_df["anomalous"],aes(anomalous)) +
#	geom_histogram(binwidth = 0.1,color ="black",fill="white")
#
#ggplot(plotting_df["vorticity"],aes(vorticity)) +
#	geom_histogram(binwidth = 1e-06,color ="black",fill="white")+
#	theme_bw()

	
#ggplot(plotting_df["surface_wind"],aes(surface_wind)) +
#	geom_histogram(binwidth = 0.25,color ="black",fill="white")+
#	theme_bw()

#ggplot(plotting_df["vor_std"],aes(vor_std)) +
#	geom_histogram(binwidth = 0.25,color ="black",fill="white")+
#	theme_bw()
#ggplot(plotting_df["T_std"],aes(T_std)) +
#	geom_histogram(binwidth = 0.25,color ="black",fill="white")+
#	theme_bw()

#ggplot(stack(plotting_df["vorticity"]), 
  #     aes(x = ind, y = values)) +geom_violin()

#library(fitdistrplus)
#T_std<-as.numeric(as.vector(unlist(plotting_df["T_std"])))
#fit.gamma <- fitdist(T_std, distr = "gamma", method = "mle")
#summary(fit.gamma)
#plot(fit.gamma)


	
	
	
#vor_std<-as.numeric(as.vector(unlist(plotting_df["vor_std"])))
#fit.normal <- fitdist(T_std, distr = "norm")
#summary(fit.normal)

#plot(fit.normal)

	
# ===============================STEPX==========================================
# Global oceanic wind speed 
# ===============================STEPX==========================================
NH_ncin<-nc_open("NH.nc")
NH<-ff(ncvar_get(NH_ncin,"ws"),dim =c(NH_ncin[["dim"]][["lon"]][["len"]],
	                  NH_ncin[["dim"]][["lat"]][["len"]],
	                  NH_ncin[["dim"]][["time"]][["len"]]),
       dimnames = list(NH_ncin[["dim"]][["lon"]][["vals"]],
                      NH_ncin[["dim"]][["lat"]][["vals"]],
                      NH_ncin[["dim"]][["time"]][["vals"]]))

SH_ncin<-nc_open("SH.nc")
SH<-ff(ncvar_get(SH_ncin,"ws"),dim =c(SH_ncin[["dim"]][["lon"]][["len"]],
	                  SH_ncin[["dim"]][["lat"]][["len"]],
	                  SH_ncin[["dim"]][["time"]][["len"]]),
       dimnames = list(SH_ncin[["dim"]][["lon"]][["vals"]],
                      SH_ncin[["dim"]][["lat"]][["vals"]],
                      SH_ncin[["dim"]][["time"]][["vals"]]))
NH_mean<-ffapply(X=NH,MARGIN = c(1,2), AFUN = mean,RETURN = TRUE)
SH_mean<-ffapply(X=SH,MARGIN = c(1,2), AFUN = mean,RETURN = TRUE)
oceanic_global_wind_speed<-mean(c(mean(NH_mean[],na.rm =TRUE),mean(SH_mean[],na.rm =TRUE)))	
	 


