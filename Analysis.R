
library(ff)
library(ffbase)
library(data.table)
library(stringr)
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
csv_files <- list.files(path = pattern="*.csv")

# Read the files in, assuming comma separator
csv_files_df <- lapply(csv_files, function(x) {
	y<-unlist(str_split(x, "[.]"))[1]
	assign(y,
	       as.ffdf(fread(x)))})

# Combine them
combined_df <- do.call("ffdfappend", lapply(csv_files_df, as.ffdf))

#combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 

nop<- which(combined_df$warm_cores[]==0 |combined_df$min_press[]==0)

combined_df <-as.ffdf(combined_df[-nop,])


nop<- which(combined_df$hourly_vertically_integrated_temp[]>= 1000)

combined_df2.0 <-as.ffdf(combined_df[-nop,])


statistic_df<- data.frame(row.names = colnames(combined_df2.0))

for(i in 1: ncol(combined_df)){
	statistic_df$mean[i]<- mean(as.numeric(unlist(combined_df2.0[i][])))
	statistic_df$std[i] <- sd(as.numeric(unlist(combined_df2.0[i][])))
	statistic_df$med[i] <- median(as.numeric(unlist(combined_df2.0[i][])))
}

combined_df2.0$vor_std<-as.ff(combined_df2.0$vorticity/statistic_df["vorticity","std"])

combined_df2.0$T_std<- as.ff(combined_df2.0[,"anomalous"]/statistic_df["anomalous","std"])

combined_df2.0$local_time<- as.ff(as.POSIXct(combined_df2.0$TIME[]*3600, origin = "1900-01-01"))

