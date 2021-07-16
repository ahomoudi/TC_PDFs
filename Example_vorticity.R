#example of vorticity
library(fields)
par(mfrow=c(3,3))
ua850<-array(runif(228*160*2,-10,10),dim = c(228,160,2))

#ua850[1:114,1:80]<- -5 ;ua850[1:114,81:160]<-5
#ua850[115:228,1:80]<- -5;ua850[115:228,81:160]<-5

image.plot(ua850[,,1])+title("u wind")

va850<-array(runif(228*160*2,-10,10),dim = c(228,160,2))


#va850[1:114,1:80]<- -5 ;va850[1:114,81:160]<- -5
#va850[115:228,1:80]<-  5;va850[115:228,81:160]<-5

image.plot(va850[,,1])+title("v wind")
longitude_v<-seq(30,228+29,1)

latitude_v<- seq(0,39.75,0.25)

DIMNAMES<- list(lon =longitude_v, lat= latitude_v,time=c(1,2))

DIM<-dim(ua850)

ws_layer_time_mean<-array(sqrt(ua850**2+va850**2),dim = DIM, dimnames = DIMNAMES)

library(reshape2)

wind.dt<-melt(ws_layer_time_mean)

colnames(wind.dt)<-c("lon","lat","time","mean_wind")

wind.dt$wind_dir<-melt(rad2deg(atan2(va850,ua850)))$value

wind.dt$u<-melt(ua850)$value
wind.dt$v<-melt(va850)$value

scaler <- 1

library(ggplot2)
p <- ggplot() +geom_raster(data=wind.dt, aes(fill=mean_wind,x=lon,y=lat),interpolate = TRUE)+
	geom_segment(data=wind.dt,
	             aes(x=lon, y=lat, xend=lon+u*scaler, yend=lat+v*scaler),
	             arrow=arrow(length = unit(0.5,"cm")), size = 0.1)+
	scale_fill_gradientn(colours = terrain.colors(5))


#print(p)

image.plot(ws_layer_time_mean[,,1])+title("wind speed")

image.plot(rad2deg(atan2(va850,ua850))[,,1])+title("wind direction")

#image(sqrt(ua850**2 + va850**2))

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

DIM<-c(dim(ua850))

#Calculation of delta x and delta y 
deltay<-vector(mode = "double",length = DIM[1])
deltax<-vector(mode = "double",length = DIM[2])

for (y in 1:DIM[2]) { 
	deltax[y]<-haversine_in_R(longitude_v[1],latitude_v[y],longitude_v[2],latitude_v[y])
}
for (x in 1:DIM[1]) { 
	deltay[x]<-haversine_in_R(longitude_v[x],latitude_v[1],longitude_v[x],latitude_v[2])
}

plot(latitude_v,deltax)+title("Delta X")

plot(longitude_v,deltay)+title("Delta Y")

deltay<-mean(deltay)

dyn.load("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/fortran_subroutines/vorticity.so")

is.loaded("relative_vorticity")

vorticity<- array(data= -99,dim = dim(va850))



#par(mfrow=c(1,3))

image.plot(vorticity[,,1])+title("empty array for vorticity")

a= 0.5

mone=DIM[1]-2 
none=DIM[2]-2
ttt<-DIM[3]
u_array<-ua850
v_array<-va850
for(tt in 1:ttt){
for (x in 2:mone){
for (y in 2:none){
	vorticity[x,y,tt] = (-1.00 *a * v_array[x,y-1,tt]/deltax[y]) 
		+ (a * v_array[x,y+1,tt]/deltax[y]) 
		- (-1.00 * a * u_array[x-1,y,tt]/deltay) 
		+ (a * u_array[x+1,y,tt]/deltay) 		
}
}
}
vorticity[vorticity==-99]<-NA

image.plot(vorticity[,,1])+title("Vorticity from R")

vorticity2<- array(data= 00,dim = dim(va850))

vorticity2<- array(.Fortran("relative_vorticity",
	    m= as.integer(DIM[1]),
	    n = as.integer(DIM[2]),
	    o= as.integer(DIM[3]),
	    u=as.numeric(ua850),
	    v=as.numeric(va850),
	    deltax =as.numeric(deltax),
	    deltay=as.numeric(deltay),
	    output_array=as.numeric(vorticity2))$output_array,
               dim = dim(ua850))

DI
image.plot(vorticity2[,,1])+title("Vorticity from fortran")

