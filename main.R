#0###############################################################################
# This code is written to calculate PDFs and statistical properties of the chosen 
#variables (vorticity, surface wind speed, and vertically integrated anomalous 
# temperature) according to :
#         Camargo, S. J., & Zebiak, S. E. (2002). Improving the Detection and 
#         Tracking of Tropical Cyclones in Atmospheric General Circulation 
#             Models. Weather and   Forecasting, 17(6), 1152–1162.
#         https://doi.org/10.1175/1520-0434(2002)017<1152:ITDATO>2.0.CO;2
#                  Written by :Ahmed Homoudi
#                      September 2020
#1.Libraries====================================================================
library(ncdf4)
library(ff)
library(ffbase)
library(reshape2)
library(tidyverse)
library(magrittr)
library(stringr)
workdir<-system("pwd",intern = TRUE)
setwd(workdir)
#2.sources======================================================================
source("functions/numextract.R")


#2.Inputs=======================================================================
# The inputs are relative vorticity at 850 hPa level, temperature at 300,500,700
#hPa in May, June, September, October and November (1990-2019)
#for the region between (0-35°N) and (30°E-80°E)

nc.files<-list.files(path = "example_data/", pattern = ".nc$")

#nc.files<-list.files(pattern = ".nc$")
array_name<-vector(length = length(nc.files))

for (i in 1:length(nc.files)) {
  ncin<- nc_open(paste0("example_data/",nc.files[i]))                           #open netcdf file 
  variable_name<-ncin[["var"]][[1]][["name"]]                                   #get variable name 
  
  year_of_data <- strsplit(unlist(str_split(nc.files[i], "_"))[5],              #get year of the data 
                           "[.]")[[1]][1]
  
  array_name[i] <- paste0(variable_name,"_",year_of_data)		#store array names in a vector 
  
  bla_bla<-numextract(unlist(str_split(nc.files[i], "_"))[3])	# extract pressure level if this is the case 
  
  if (!is.na(bla_bla)){

      if(any(bla_bla==c("300","500","700","850"))){
        
        press_lev  <- numextract(unlist(str_split(nc.files[i], "_"))[3])
        array_name[i] <- paste0(variable_name,press_lev,"_",year_of_data)	# modifiy array name in case of press level
      }
  }
    

  med<-ncvar_get(ncin,variable_name)			# get the variable 
  longitude_v<-ncvar_get(ncin,"longitude")
  latitude_v<-ncvar_get(ncin,"latitude")
  raw_time<-ncvar_get(ncin,"time")
  
  assign(array_name[i],ff(med,dimnames =  
                         list(longitude_v,latitude_v,raw_time),dim = dim(med))) #save variable array on the disk
  
  rm(med,ncin)				# remove dummy variables 
  
}
# ==============================STEP1===========================================
#Hourly vertically integrated temperature, the sum of the temperature at 
#         700, 500, and 300 hPa in each grid point.
temperature_on_levels <- grep(pattern = '00', x = array_name,value = TRUE)	#get temp from array name vector 

hourly_vertically_integrated_temp<-ff(array(data = 0.0),dim = 
	                  	dim(get(temperature_on_levels[1])))	#dummy variable 

for (i in  1:length(temperature_on_levels)){
  hourly_vertically_integrated_temp<-ff(hourly_vertically_integrated_temp[]+
    get(temperature_on_levels[i])[],
    dim = dim(get(temperature_on_levels[1])),
    dimnames = dimnames(get(temperature_on_levels[1])))		#cumulative sum of temperature 
}
# ==============================STEP2===========================================
# A local average hourly vertically integrated temperature is calculated in a 
#   square, 7X7 grid points, centred at the grid point of interest.

dyn.load("fortran_subroutines/vert_int_temp.so")		#load Fortran subroutine 

is.loaded("integrated_temperature")			#check if it is loaded 

hourly_vertically_integrated_temp_mean<-ff(.Fortran("integrated_temperature",
                                                       input_array =as.numeric(hourly_vertically_integrated_temp[]),
                                                       lon = as.integer(length(longitude_v)),
                                                       lat = as.integer(length(latitude_v)),
                                                       time = as.integer(length(raw_time)),
                                                       output_array=as.numeric(hourly_vertically_integrated_temp[]))$output_array,
                                              dim = dim(hourly_vertically_integrated_temp),
                                              dimnames =dimnames(hourly_vertically_integrated_temp))

hourly_vertically_integrated_temp_df<-reshape2::melt(hourly_vertically_integrated_temp[],
                                                     varnames = names(dimnames(hourly_vertically_integrated_temp)),
                                                     value.name="hourly_vertically_integrated_temp")%>%as.ffdf()

colnames(hourly_vertically_integrated_temp_df)<-c("LON","LAT","TIME","hourly_vertically_integrated_temp")


hourly_vertically_integrated_temp_mean_df<-reshape2::melt(hourly_vertically_integrated_temp_mean[],
                                                     varnames = names(dimnames(hourly_vertically_integrated_temp_mean)),
                                                     value.name="hourly_vertically_integrated_temp_mean")%>%as.ffdf()

colnames(hourly_vertically_integrated_temp_mean_df)<-c("LON","LAT","TIME","hourly_vertically_integrated_temp_mean")


hourly_vertically_integrated_temp_df$hourly_vertically_integrated_temp_mean<-
  hourly_vertically_integrated_temp_mean_df$hourly_vertically_integrated_temp_mean



# ==============================STEP3===========================================
# The anomalous vertically integrated hourly temperature in the grid point at the 
#center of the square is calculated as the difference between the 
#hourly vertically integrated temperature in that grid point and the 
#local average hourly vertically integrated temperature around it.

hourly_vertically_integrated_temp_df$anomalous <- hourly_vertically_integrated_temp_df$hourly_vertically_integrated_temp -
  hourly_vertically_integrated_temp_df$hourly_vertically_integrated_temp_mean

rm(hourly_vertically_integrated_temp,
   hourly_vertically_integrated_temp_mean,
   hourly_vertically_integrated_temp_mean_df)

# ==============================STEP4===========================================  
#The anomalous hourly temperature for the levels 850, 700, 500, and 300 hPa is 
#calculated, by analogy to the anomalous vertically integrated temperature.

#AT 300 hPa#####################################################################
temp300<-grep(pattern = '300', x = array_name,value = TRUE)

dyn.load("fortran_subroutines/vert_int_temp.so")

is.loaded("integrated_temperature")

hourly_square_temp_mean<-ff(.Fortran("integrated_temperature",
                                                    input_array =as.numeric(get(temp300)[]),
                                                    lon = as.integer(length(longitude_v)),
                                                    lat = as.integer(length(latitude_v)),
                                                    time = as.integer(length(raw_time)),
                                                    output_array=as.numeric(get(temp300)[]))$output_array,
                                           dim = dim(get(temp300)[]),
                                           dimnames =dimnames(get(temp300)))

hourly_square_temp_mean_df<-reshape2::melt(hourly_square_temp_mean[],
                                                     varnames = names(dimnames(hourly_square_temp_mean)),
                                                     value.name="hourly_temp300_mean")%>%as.ffdf()

colnames(hourly_square_temp_mean_df)<-c("LON","LAT","TIME","temp300_mean")

hourly_vertically_integrated_temp_df$temp300_mean<-hourly_square_temp_mean_df$temp300_mean

hourly_temp_df<-reshape2::melt(get(temp300)[],
                               varnames = names(dimnames(get(temp300))),
                               value.name="temp300")%>%as.ffdf()

hourly_vertically_integrated_temp_df$temp300<-hourly_temp_df$temp300

hourly_vertically_integrated_temp_df$temp300_anomaly<-hourly_vertically_integrated_temp_df$temp300-
  hourly_vertically_integrated_temp_df$temp300_mean
rm(hourly_square_temp_mean,
   hourly_temp_df,
   hourly_square_temp_mean_df)

#AT 500 hPa#####################################################################
temp500<-grep(pattern = '500', x = array_name,value = TRUE)

dyn.load("fortran_subroutines/vert_int_temp.so")

is.loaded("integrated_temperature")

hourly_square_temp_mean<-ff(.Fortran("integrated_temperature",
                                     input_array =as.numeric(get(temp500)[]),
                                     lon = as.integer(length(longitude_v)),
                                     lat = as.integer(length(latitude_v)),
                                     time = as.integer(length(raw_time)),
                                     output_array=as.numeric(get(temp500)[]))$output_array,
                            dim = dim(get(temp500)[]),
                            dimnames =dimnames(get(temp500)))

hourly_square_temp_mean_df<-reshape2::melt(hourly_square_temp_mean[],
                                           varnames = names(dimnames(hourly_square_temp_mean)),
                                           value.name="hourly_temp500_mean")%>%as.ffdf()

colnames(hourly_square_temp_mean_df)<-c("LON","LAT","TIME","temp500_mean")

hourly_vertically_integrated_temp_df$temp500_mean<-hourly_square_temp_mean_df$temp500_mean

hourly_temp_df<-reshape2::melt(get(temp500)[],
                               varnames = names(dimnames(get(temp500))),
                               value.name="temp500")%>%as.ffdf()

hourly_vertically_integrated_temp_df$temp500<-hourly_temp_df$temp500

hourly_vertically_integrated_temp_df$temp500_anomaly<-hourly_vertically_integrated_temp_df$temp500-
  hourly_vertically_integrated_temp_df$temp500_mean

rm(hourly_square_temp_mean,
   hourly_temp_df,
   hourly_square_temp_mean_df)
#AT 700 hPa#####################################################################
temp700<-grep(pattern = '700', x = array_name,value = TRUE)

dyn.load("fortran_subroutines/vert_int_temp.so")

is.loaded("integrated_temperature")

hourly_square_temp_mean<-ff(.Fortran("integrated_temperature",
                                     input_array =as.numeric(get(temp700)[]),
                                     lon = as.integer(length(longitude_v)),
                                     lat = as.integer(length(latitude_v)),
                                     time = as.integer(length(raw_time)),
                                     output_array=as.numeric(get(temp700)[]))$output_array,
                            dim = dim(get(temp700)[]),
                            dimnames =dimnames(get(temp700)))

hourly_square_temp_mean_df<-reshape2::melt(hourly_square_temp_mean[],
                                           varnames = names(dimnames(hourly_square_temp_mean)),
                                           value.name="hourly_temp700_mean")%>%as.ffdf()

colnames(hourly_square_temp_mean_df)<-c("LON","LAT","TIME","temp700_mean")

hourly_vertically_integrated_temp_df$temp700_mean<-hourly_square_temp_mean_df$temp700_mean

hourly_temp_df<-reshape2::melt(get(temp700)[],
                               varnames = names(dimnames(get(temp700))),
                               value.name="temp700")%>%as.ffdf()

hourly_vertically_integrated_temp_df$temp700<-hourly_temp_df$temp700

hourly_vertically_integrated_temp_df$temp700_anomaly<-hourly_vertically_integrated_temp_df$temp700-
  hourly_vertically_integrated_temp_df$temp700_mean

rm(hourly_square_temp_mean,
   hourly_temp_df,
   hourly_square_temp_mean_df)

#AT 850 hPa#####################################################################
temp850<-grep(pattern = '850', x = array_name,value = TRUE)

dyn.load("fortran_subroutines/vert_int_temp.so")

is.loaded("integrated_temperature")

hourly_square_temp_mean<-ff(.Fortran("integrated_temperature",
                                     input_array =as.numeric(get(temp850)[]),
                                     lon = as.integer(length(longitude_v)),
                                     lat = as.integer(length(latitude_v)),
                                     time = as.integer(length(raw_time)),
                                     output_array=as.numeric(get(temp850)[]))$output_array,
                            dim = dim(get(temp850)[]),
                            dimnames =dimnames(get(temp850)))

hourly_square_temp_mean_df<-reshape2::melt(hourly_square_temp_mean[],
                                           varnames = names(dimnames(hourly_square_temp_mean)),
                                           value.name="hourly_temp850_mean")%>%as.ffdf()

colnames(hourly_square_temp_mean_df)<-c("LON","LAT","TIME","temp850_mean")

hourly_vertically_integrated_temp_df$temp850_mean<-hourly_square_temp_mean_df$temp850_mean

hourly_temp_df<-reshape2::melt(get(temp850)[],
                               varnames = names(dimnames(get(temp850))),
                               value.name="temp850")%>%as.ffdf()

hourly_vertically_integrated_temp_df$temp850<-hourly_temp_df$temp850

hourly_vertically_integrated_temp_df$temp850_anomaly<-hourly_vertically_integrated_temp_df$temp850-
  hourly_vertically_integrated_temp_df$temp850_mean

rm(hourly_square_temp_mean,
   hourly_temp_df,
   hourly_square_temp_mean_df)
# ==============================STEP5&6===========================================  
#The difference of the temperature anomalies at 850 and 300 hPa is calculated. 
#The temperature anomaly at 850 hPa has to be positive and smaller than the 
#temperature anomaly at 300 hPa. The signs of the temperature anomalies at 700,
#500, and 300 hPa are compared.

dyn.load("fortran_subroutines/warm_core_subroutine.so")

is.loaded("warm_core")

warm_core_filter<- ff(array(data= 0.00),dim = dim(get(temp300)))
warm_core_filter<-ff(.Fortran("warm_core",
                              lon = as.integer(length(longitude_v)),
                              lat = as.integer(length(latitude_v)),
                              time = as.integer(length(raw_time)),
                              temp300=as.numeric(array(hourly_vertically_integrated_temp_df$temp300_anomaly,dim =dim(get(temp300)[]))[]),
                              temp500=as.numeric(array(hourly_vertically_integrated_temp_df$temp500_anomaly,dim =dim(get(temp500)[]))[]),
                              temp700=as.numeric(array(hourly_vertically_integrated_temp_df$temp700_anomaly,dim =dim(get(temp700)[]))[]),
                              temp850=as.numeric(array(hourly_vertically_integrated_temp_df$temp850_anomaly,dim =dim(get(temp850)[]))[]),
                              filter=as.numeric(warm_core_filter[]))$filter,
                     dim = dim(get(temp850)[]),
                     dimnames =dimnames(get(temp850)))


warm_core_df<-reshape2::melt(warm_core_filter[],value.name="warm_core_filter")%>%as.ffdf()

colnames(warm_core_df)<-c("LON","LAT","TIME","warm_cores")

hourly_vertically_integrated_temp_df$warm_cores<-warm_core_df$warm_cores

rm(warm_core_df,warm_core_filter)

temperature_on_levels <- grep(pattern = 't', x = array_name,value = TRUE)	#get temp from array name vector

rm(list = temperature_on_levels)

# ==============================STEP7===========================================  
#The sea level pressure is the minimum in a centred 7X7 box.
surface_pressure<-grep(pattern = 'sp', x = array_name,value = TRUE)

dyn.load("fortran_subroutines/pressure_minima.so")

is.loaded("minimum_press")

min_press_filter<- ff(array(data= 0.00),dim = dim(get(surface_pressure)))
min_press_filter<-ff(.Fortran("minimum_press",
                              lon = as.integer(length(longitude_v)),
                              lat = as.integer(length(latitude_v)),
                              time = as.integer(length(raw_time)),
                              a=as.integer(7),
                              b=as.integer(7),
                              ps=as.numeric(get(surface_pressure)[]),
                              box= as.numeric(array(0,dim = c(7,7))),
                              filter=as.numeric(min_press_filter[]))$filter,
                     dim = dim(get(surface_pressure)),
                     dimnames =dimnames(get(surface_pressure)))

min_press_df<-reshape2::melt(min_press_filter[],value.name="min_press_filter")%>%as.ffdf()

colnames(min_press_df)<-c("LON","LAT","TIME","min_press")

hourly_vertically_integrated_temp_df$min_press<-min_press_df$min_press

rm(min_press_df,min_press_filter)

rm(list= surface_pressure)
# ==============================STEP8===========================================  
# Obtain surface wind 
u10<-grep(pattern = 'u10', x = array_name,value = TRUE)
v10<-grep(pattern = 'v10', x = array_name,value = TRUE)

ws<-ff(sqrt((get(u10)[])**2 +(get(v10)[])**2),
      dim = dim(get(u10)),
      dimnames = dimnames(get(u10)))

surface_wind_df<-reshape2::melt(ws[],value.name="surface_wind")%>%as.ffdf()

colnames(surface_wind_df)<-c("LON","LAT","TIME","surface_wind")

hourly_vertically_integrated_temp_df$surface_wind<-surface_wind_df$surface_wind

rm(ws,surface_wind_df,list = u10)
rm(list = v10)

#Obtain vorticity 
vorticity <-grep(pattern = 'vo', x = array_name,value = TRUE)


vorticity_df<-reshape2::melt(get(vorticity)[],value.name="vorticity")%>%as.ffdf()

colnames(vorticity_df)<-c("LON","LAT","TIME","vorticity")

hourly_vertically_integrated_temp_df$vorticity<-vorticity_df$vorticity

rm(vorticity_df)

rm(list = vorticity)
# ==============================STEP9===========================================  
#Filtering and writing to text files 
nop <- which(hourly_vertically_integrated_temp_df$warm_core[]==0.0 &
               hourly_vertically_integrated_temp_df$min_press[] == 0.0)

hourly_vertically_integrated_temp_df_filtered<- hourly_vertically_integrated_temp_df[-nop,]%>%as.ffdf()


print(c("The number of unfilttered points is ", nrow(hourly_vertically_integrated_temp_df)))

print(c("The number of filttered points is ", nrow(hourly_vertically_integrated_temp_df_filtered)))

write.csv(hourly_vertically_integrated_temp_df_filtered[],file =paste0("warm_core_points_",year_of_data,"_.txt"))

write.csv(hourly_vertically_integrated_temp_df[],file =paste0("all_grid_with_filtervalues_",year_of_data,"_.txt"))



