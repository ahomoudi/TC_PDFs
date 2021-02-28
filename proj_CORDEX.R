# function  to read opened  COREDEX netcdf file and projected to global CRS + crop it 

proj_read_toff<- function(netcdf_file){
	library(ncdf4)
	library(stringr)
	library(ff)
	
	ncin <- nc_open(netcdf_file)
	
	file_name<-ncin[["filename"]]
	
	variable_name<- paste0(ncin[["var"]][[2]][["name"]],"_",
		   strsplit(unlist(str_split(file_name,
		   	      "_"))[9],     #get year of the data 
		            "[.]")[[1]][1])

	nc.var<- ncin[["var"]][[2]][["name"]]		#varaible name 
	
	
	output_file<- paste0(unlist(strsplit(file_name,"[.]"))[1],  #projected file name
		 "_tmp.nc")
	
	first_command<-paste0("cdo remapbil,proj_WAS.txt ",        #cdo command 
		  file_name," ",
		  output_file)
	
	system(first_command)		# call cdo in the system 
	
	ncin_cropped<- nc_open(output_file)		#open netcdf file 
	
	med<- ncvar_get(ncin_cropped,nc.var)		#obtain variable from netcdf file 
	
	lon <-ncvar_get(ncin_cropped,"longitude")
	lat <- ncvar_get(ncin_cropped,"latitude")
	time <- ncvar_get(ncin_cropped,"time")
	
	processed_array<- ff(med, dim= dim(med),
		 dimnames= list(lon,lat,time))	# save data as ff object
	
	system(paste0("rm ",output_file))		#remove projected netcdf file 
	
	return(processed_array)
	
}
#netcdf_file<-nc.files[i]


