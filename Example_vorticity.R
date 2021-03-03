#example of vorticity

u850<-get(ua850)[1:4,1:4,1]

v850<-get(va850)[1:4,1:4,1]

lon = lon_lonlat[1:40,1:40]

lat = lat_lonlat[1:40,1:40]

image(sqrt(u850**2 + v850**2))

haversine_in_R<-function(lon1,lat1,lon2,lat2){
	library(REdaS)
	R_EARTH<- 6371000  # radius of Earth in meters
	phi_1 <- deg2rad(lat1)
	phi_2 <- deg2rad(lat2)
	
	delta_phi <- deg2rad(lat2-lat1)
	delta_lambda <- deg2rad(lon2-lon1)
	
	a <- sin(delta_phi / 2.0) ** 2 + cos(phi_1) * cos(phi_2) * sin(delta_lambda / 2.0) ** 2
	
	c = 2 * atan2(sqrt(a),sqrt(1 - a))
	
	meters = R_EARTH * c  # output distance in meters
	return(meters)
}

delta<-array(00, dim = c(365,239,2))

for ( i in 2:240){ # in lat
	
	lon_col<-lon_lonlat[,i-1]
	lat_col<-lat_lonlat[,i-1]
	result_col<-lat_col
	for ( j in 2:366){
	result_col[j-1]<-haversine_in_R(lon1 = lon_col[j-1] ,lon2 = lon_col[j] ,
		            lat1 =lat_col[j-1], lat2 =lat_col[j])	
	}
	delta[,i-1,1]<-result_col[-length(result_col)]
 
}

for ( i in 2:366){ # in lat
	
	lon_row<-lon_lonlat[i-1,]
	lat_row<-lat_lonlat[i-1,]
	result_row<-lat_row
	for ( j in 2:240){
		result_row[j-1]<-haversine_in_R(lon1 = lon_col[j-1] ,lon2 = lon_col[j] ,
			            lat1 =lat_col[j-1], lat2 =lat_col[j])	
	}
	delta[i-1,,2]<-result_row[-length(result_row)]
	
}

library(fields)
par(mfrow=c(1,2))

image.plot(delta[,,1])

image.plot(delta[,,2])

image.plot(delta[,,1]*delta[,,2])
