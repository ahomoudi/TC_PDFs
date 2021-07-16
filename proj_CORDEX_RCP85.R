# function  to read opened  COREDEX netcdf file and projected to global CRS + crop it 

proj_read_toff<- function(netcdf_file,refrence_file){
	#library(ncdf4)
	library(stringr)
	#library(ff)
	#library(berryFunctions)
	
	ncin <- nc_open(netcdf_file)
	
	file_name<-ncin[["filename"]]
	
	variable_name<- paste0(ncin[["var"]][[2]][["name"]],"_",
		   strsplit(unlist(str_split(file_name,
		   	      "_"))[9],     #get year of the data 
		   	      "[.]")[[1]][1])
	
	nc.var<- ncin[["var"]][[2]][["name"]]		#varaible name 
	
	
	output_file<- paste0(unlist(strsplit(file_name,"[.]"))[1],  #projected file name
		 "_proj.nc")
	
	
	if(all(nc.var!=c("ua850","va850","uas","vas","ua300","va300"))){
	
	    first_command<-paste0("cdo remapbil,proj_WAS.txt ",        #cdo command 
			  file_name," ",
			  output_file)
	    
	    system(first_command)		# re-projection 
	    
	    system(paste0("mv ",output_file," projected/",
	                  output_file))		#remove projected netcdf file 
	}

	# this part with the help from Dr. Silje Sorland from ETH Zurich
	if(any(nc.var==c("ua850","va850","ua300","va300"))){		# for wind over pressure levels
		
	   output_file1<- paste0(unlist(strsplit(file_name,"[.]"))[1],  #tmp file name
			 "_tmp.nc")
		
	   first_command<-paste0("ncks -C -v ",nc.var,",time,rotated_pole,plev ",file_name," ",output_file1)
	   
	   #Only copying over the variables that are correct to a new file
	   system(first_command)
	   
	   second_command<-paste0("ncks -A -v rlat ",refrence_file," ",output_file1)
	   
	   system(second_command)
	   
	   third_command<-paste0("ncks -A -v rlon ",refrence_file," ",output_file1)
	   
	   system(third_command)
	   
	   fourth_command<-paste0("ncks -A -v lon ",refrence_file," ",output_file1)
	   
	   system(fourth_command)
	   
	   fifth_command<-paste0("ncks -A -v lat ",refrence_file," ",output_file1)
	   
	   system(fifth_command)
	   
	   sixth_command<-paste0("cdo remapbil,proj_WAS.txt ",        #cdo command 
	   	  output_file1," ",
	   	  output_file)
	   
	   system(sixth_command)		# re-projection
	   
	   system(paste0("mv ",output_file," projected/",
	                 output_file))
	   
	   system(paste0("rm ",output_file1))
	   
	   
	}

	if(any(nc.var==c("uas","vas"))){		# for surface wind 

				
 	 output_file1<- paste0(unlist(strsplit(file_name,"[.]"))[1],  #projected file name
					  "_tmp.nc")
				
	 first_command<-paste0("ncks -C -v ",nc.var,",time,rotated_pole,height ",file_name," ",output_file1)
				
	#Only copying over the variables that are correct to a new file
	  system(first_command)
				
	second_command<-paste0("ncks -A -v rlat ",refrence_file," ",output_file1)
				
	system(second_command)
				
	third_command<-paste0("ncks -A -v rlon ",refrence_file," ",output_file1)
				
	system(third_command)
				
	fourth_command<-paste0("ncks -A -v lon ",refrence_file," ",output_file1)
				
	system(fourth_command)
				
	fifth_command<-paste0("ncks -A -v lat ",refrence_file," ",output_file1)
				
	system(fifth_command)
	
	#cdo command 			
	sixth_command<-paste0("cdo remapbil,proj_WAS.txt ",output_file1," ",output_file)
				
	system(sixth_command)		# re-projection 
				
				
	system(paste0("mv ",output_file," projected/",output_file))	   
	
	system(paste0("rm ",output_file1))
				
	}

}

