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
library(REdaS)
#workdir<-system("pwd",intern = TRUE)
#setwd(workdir)

options(fftempdir = "/media/ahmed/Daten/ff")

getOption("fftempdir")

#2.sources======================================================================
source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/functions/numextract.R")
source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/functions/haversine.R")
source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/proj_CORDEX.R")

setwd("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/cordex")
#2.Inputs=======================================================================
# The inputs are relative vorticity at 850 hPa level, temperature at 300,500,700
#hPa in May, June, September, October and November (1990-2019)
#for the region between (0-35°N) and (30°E-80°E)

nc.files<-list.files(pattern = ".nc$")
nc.files<-nc.files[-c(1)]

nc.files

array_name<-vector(length = length(nc.files))

for (i in 1:length(nc.files)) {
  
  ncin<- nc_open(nc.files[i])                                                  #open netcdf file 
  #ncin<- nc_open(paste0("example_data/",nc.files[i])) 
  variable_name<-ncin[["var"]][[2]][["name"]]                                   #get variable name 
  # 1 for ERA5 data 
  year_of_data <- strsplit(unlist(str_split(nc.files[i], "_"))[9],              #get year of the data 
                           "[.]")[[1]][1]
  #6 instead of 5 for 1990/6 for ERA5 1991-2019
  
  array_name[i] <- paste0(variable_name,"_",year_of_data)		#store array names in a vector 
}

rm(ncin)				          # remove dummy variables 
# CORDEX data  ------------------------------------------------------------
system("mkdir projected")

for (i in 1:length(nc.files)) {
  
  ta850 <-grep(pattern = 'ta850', x = nc.files,value = TRUE)
  
  print(i)
  assign(array_name[i],proj_read_toff(netcdf_file=nc.files[i],
                                      refrence_file= ta850))
  
  #rm(ncin)				          # remove dummy variables 
  #rm(med)
  gc()
  
}

# MSL calculation  ----------------------------------------------------------

#The calculate mean sea level pressure from geopotinal height and temperature at 850 hpa 
geopotanisal_height <-grep(pattern = 'zg', x = array_name,value = TRUE)

temp850<-grep(pattern = 'ta850', x = array_name,value = TRUE)

dim(get(geopotanisal_height))
dim(get(temp850))

DIMnames<-dimnames(get(temp850))

longitude_v<-DIMnames[[1]]
latitude_v<-DIMnames[[2]]
raw_time<-DIMnames[[3]]

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/msl_calc.so")

is.loaded("p850_to_msl")

msl_array<- ff(array(data= 0.00),dim = dim(get(geopotanisal_height)))

msl_array<-ff(.Fortran("p850_to_msl",
                       lon = as.integer(length(longitude_v)),
                       lat = as.integer(length(latitude_v)),
                       time = as.integer(length(raw_time)),
                       zg_array=as.numeric(get(geopotanisal_height)[]),
                       ta_array=as.numeric(get(temp850)[]),
                       output_array=as.numeric(msl_array[]))$output_array,
              dim = dim(get(geopotanisal_height)),
              dimnames =dimnames(get(geopotanisal_height)))

#calculate vorticity --------------------------------------------------------
ua850 <-grep(pattern = 'ua850', x = array_name,value = TRUE)

va850 <-grep(pattern = 'va850', x = array_name,value = TRUE)

DIMnames<-dimnames(get(ua850))

longitude_v<-DIMnames[[1]]
latitude_v<-DIMnames[[2]]
raw_time<-DIMnames[[3]]

DIM<-c(length(longitude_v),length(latitude_v),length(raw_time))

#Calculation of delta x and delta y 
deltay<-vector(mode = "double",length = DIM[1])
deltax<-vector(mode = "double",length = DIM[2])

for (y in 1:DIM[2]) { 
  deltax[y]<-haversine_in_R(longitude_v[50],latitude_v[y],longitude_v[51],latitude_v[y])
}
for (x in 1:DIM[1]) { 
  deltay[x]<-haversine_in_R(longitude_v[x],latitude_v[35],longitude_v[x],latitude_v[36])
}

deltay<-mean(deltay)

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vorticity.so")

is.loaded("relative_vorticity")


vorticity<- ff(array(data= -99),dim = dim(get(va850)))

vorticity<- ff(.Fortran("relative_vorticity",
                        m= as.integer(DIM[1]),
                        n = as.integer(DIM[2]),
                        o= as.integer(DIM[3]),
                        u=as.numeric(get(ua850)[]),
                        v=as.numeric(get(va850)[]),
                        deltax =as.numeric(deltax),
                        deltay=as.numeric(deltay),
                        output_array=as.numeric(vorticity[]))$output_array,
               dim = dim(get(ua850)),
               dimnames =dimnames(get(ua850)))

#vorticity[vorticity==-99]<-NA

rm(ua850,va850)
# ==============================STEP1===========================================
#Hourly vertically integrated temperature, the sum of the temperature at 
#         700, 500, and 300 hPa in each grid point.

temperature_on_levels<-c(grep(pattern = '300', x = array_name,value = TRUE),
                         grep(pattern = '500', x = array_name,value = TRUE),
                         grep(pattern = '700', x = array_name,value = TRUE))

#temperature_on_levels <- grep(pattern = '00', x = array_name,value = TRUE)	#get temp from array name vector 

hourly_vertically_integrated_temp<-ff(array(data = 0.0),dim = 
	                  	dim(get(temperature_on_levels[1])))	#dummy variable 

for (i in  1:length(temperature_on_levels)){
  hourly_vertically_integrated_temp<-ff(hourly_vertically_integrated_temp[]+
    get(temperature_on_levels[i])[],
    dim = dim(get(temperature_on_levels[1])),
    dimnames = dimnames(get(temperature_on_levels[1])))		#cumulative sum of temperature 
}
#rm()
# ==============================STEP2===========================================
# A local average hourly vertically integrated temperature is calculated in a 
#   square, 7X7 grid points, centred at the grid point of interest.

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vert_int_temp.so")		#load Fortran subroutine 

is.loaded("integrated_temperature")			#check if it is loaded 

DIMnames<-dimnames(get(temperature_on_levels[2]))

longitude_v<-DIMnames[[1]]
latitude_v<-DIMnames[[2]]
raw_time<-DIMnames[[3]]


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

#gc()

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

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vert_int_temp.so")

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

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vert_int_temp.so")

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

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vert_int_temp.so")

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

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vert_int_temp.so")

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

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/warm_core_subroutine.so")

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

# find the eye of the storm ----------------------------------------------
dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/pressure_minimum/lowest_pressure_eye_storm.so")

is.loaded("eyeofthestorm")

min_press_filter<- ff(array(data= 0.00),dim = dim(get(geopotanisal_height)))

min_press_filter<-ff(.Fortran("eyeofthestorm",
                              lon = as.integer(length(longitude_v)),
                              lat = as.integer(length(latitude_v)),
                              time = as.integer(length(raw_time)),
                              a=as.integer(7),
                              b=as.integer(7),
                              ps=as.numeric(msl_array[]),
                              box= as.numeric(array(0,dim = c(7,7))),
                              filter=as.numeric(min_press_filter[]))$filter,
                     dim = dim(get(geopotanisal_height)),
                     dimnames =dimnames(get(geopotanisal_height)))

gc()

min_press_df<-reshape2::melt(min_press_filter[],value.name="min_press_filter")%>%as.ffdf()

colnames(min_press_df)<-c("LON","LAT","TIME","min_press")

hourly_vertically_integrated_temp_df$min_press<-min_press_df$min_press

rm(min_press_df,min_press_filter)

#rm(list= surface_pressure)
# ==============================STEP8===========================================  
# Obtain surface wind 
u10<-grep(pattern = 'uas', x = array_name,value = TRUE)
v10<-grep(pattern = 'vas', x = array_name,value = TRUE)

ws<-ff(sqrt((get(u10)[])**2 +(get(v10)[])**2),
      dim = dim(get(u10)),
      dimnames = dimnames(get(u10)))

surface_wind_df<-reshape2::melt(ws[],value.name="surface_wind")%>%as.ffdf()

colnames(surface_wind_df)<-c("LON","LAT","TIME","surface_wind")

hourly_vertically_integrated_temp_df$surface_wind<-surface_wind_df$surface_wind

rm(ws,surface_wind_df,list = u10)
rm(list = v10)

# Obtain Vorticity --------------------------------------------------------

vorticity_df<-reshape2::melt(vorticity[],value.name="vorticity")%>%as.ffdf()

colnames(vorticity_df)<-c("LON","LAT","TIME","vorticity")

hourly_vertically_integrated_temp_df$vorticity<-vorticity_df$vorticity

rm(vorticity_df)

# ==============================STEP9===========================================  
#Filtering and writing to text files 
nop <- which(hourly_vertically_integrated_temp_df$warm_core[]==0.0 &
               hourly_vertically_integrated_temp_df$min_press[] == 0.0)

hourly_vertically_integrated_temp_df_filtered<- hourly_vertically_integrated_temp_df[-nop,]%>%as.ffdf()


print(c("The number of unfilttered points is ", nrow(hourly_vertically_integrated_temp_df)))

print(c("The number of filttered points is ", nrow(hourly_vertically_integrated_temp_df_filtered)))


#write.table.ffdf(x= hourly_vertically_integrated_temp_df_filtered,
               #  file =paste0("/lustre/scratch2/ws/1/ahho623a-FRM_PDFs_project/Results/warm_core_points_",year_of_data,"_.txt"))

write.table.ffdf(x= hourly_vertically_integrated_temp_df_filtered,
                 file =paste0("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/Results/warm_core_points_",year_of_data,"_.txt"))


#write.table.ffdf(x= hourly_vertically_integrated_temp_df,
                 #file =paste0("/lustre/scratch2/ws/1/ahho623a-FRM_PDFs_project/Results/all_grid_with_filtervalues_",year_of_data,"_.txt"))



# ==============================END===========================================  