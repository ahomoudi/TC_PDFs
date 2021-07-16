
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


points_processed <- read.csv("/media/ahmed/Volume/TC_FRM/IBTrACS.NI.list.v04r00.lines/points_processed.csv")


# Read function -----------------------------------------------------------

source("/media/ahmed/Daten/GCMs/ITCZ_analysis/extra_function.R")
source("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/points_to_line.R")



# Process -----------------------------------------------------------------

## example data
d <- points_processed

coordinates(d) <- ~LON+LAT

## list of Lines per id, each with one Line in a list
x <- lapply(split(d, d$Storm), function(x) Lines(list(Line(coordinates(x))), x$Storm[1L]))

# the corrected part goes here:
lines <- SpatialLines(x)
data <- data.frame(id = unique(d$Storm))
rownames(data) <- data$id
l <- SpatialLinesDataFrame(lines, data)

proj4string(l) <- CRS( "+proj=longlat +datum=WGS84" )

# Write to shape file 
writeOGR(obj=l, dsn="tempdir", layer="lines", driver="ESRI Shapefile") 

# plot tracks 
#leaflet(data=lines) %>% addTiles() %>% addPolylines()

mapview(l)

l_fortify <- fortify(l)


ggplot(l_fortify, aes(x=long, y=lat, group=group)) + 
	geom_path() +
	theme_classic()


# Features analysis ------------------------------------------------------------

# for vorticity----------------------------------------------------------------

Results<-data.frame(NO=seq(1,30,1))

BREAKS<-seq(from = 0.00, to = 0.003, by = 1e-04)

Results$Vorticity<-seq(from = 5e-05, to = 0.00295, by = 1e-04)

Results$Vorticity_PDF_all<-table(cut(points_processed$vo850,BREAKS, right=FALSE)) / length(points_processed$vo850) 

Vorticity_PDF_MTS<-points_processed$vo850[-which(points_processed$MTS=='no')]

Results$Vorticity_PDF_MTS<-table(cut(Vorticity_PDF_MTS,BREAKS, right=FALSE)) / length(Vorticity_PDF_MTS) 


Fig_1a<-ggplot(Results,aes(x =Vorticity*1e03))+
	geom_path( aes(y = Vorticity_PDF_all,color = "All"))+
	geom_line( aes(y = Vorticity_PDF_MTS,color = "MTS"))+

	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(-2.5,2.5), expand = c(0, 0),
	                   breaks = seq(from = -2.5, to = 2.5, by = 0.25),
	                   labels = every_nth(seq(from = -2.5, to = 2.5, by = 0.25),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.3), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.3, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.3, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = 'bottom',
	      legend.title = element_blank(),
	      legend.text = element_text(size = 14),
	      legend.box.background = element_blank(),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Relative Vorticity "(s^{-1}~x10^{3})))+
	ylab("Relative Frequency")+
	labs(title = expression("(a)"))+
	theme(axis.text = element_text(colour = "black"))

# for eye pressure -------------------------------------------------------------

Results1<-data.frame(NO=seq(1,140,1))

BREAKS<-seq(from = 950*100, to = 1020*100, by = 50)

Results1$Pressure<-seq(from = 95025, to = 101975, by = 50)

Results1$Pressure_PDF_all<-table(cut(points_processed$Eye_Press,BREAKS, right=FALSE)) / length(points_processed$Eye_Press) 

Pressure_PDF_MTS<-points_processed$Eye_Press[-which(points_processed$MTS=='no')]

Results1$Pressure_PDF_MTS<-table(cut(Pressure_PDF_MTS,BREAKS, right=FALSE)) / length(Pressure_PDF_MTS) 


Fig_1b<-ggplot(Results1,aes(x =Pressure/100))+
	geom_path( aes(y = Pressure_PDF_all,color = "All"))+
	geom_line( aes(y = Pressure_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(950,1020), expand = c(0, 0),
	                   breaks = seq(from = 950, to = 1020, by = 5),
	                   labels = every_nth(seq(from = 950, to = 1020, by = 5),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.05), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.05, by = 0.005),
	                   labels = every_nth(seq(from = 0, to = 0.05, by = 0.005),2, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Eye Mean Sea Level Pressure [hPa]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(b)"))+
	theme(axis.text = element_text(colour = "black"))

# for wind speed   -------------------------------------------------------------

BREAKS<-seq(from = 0.00, to = 30, by = 1)

Results$ws10<-seq(from = 0.5, to = 29.5, by = 1)

Results$ws10_PDF_all<-table(cut(points_processed$w10,BREAKS, right=FALSE)) / length(points_processed$w10) 

ws10_PDF_MTS<-points_processed$w10[-which(points_processed$MTS=='no')]

Results$ws10_PDF_MTS<-table(cut(ws10_PDF_MTS,BREAKS, right=FALSE)) / length(ws10_PDF_MTS) 


Fig_1c<-ggplot(Results,aes(x =ws10))+
	geom_path( aes(y = ws10_PDF_all,color = "All"))+
	geom_line( aes(y = ws10_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,40), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 40, by = 5),
	                   labels = every_nth(seq(from = 0, to = 40, by = 5),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.3), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.3, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.3, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Wind Speed [m/s]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(c)"))+
	theme(axis.text = element_text(colour = "black"))

# for maximum sustained  wind speed   -------------------------------------------------------------

Results2<-data.frame(NO=seq(1,40,1))

BREAKS<-seq(from = 0.00, to = 40, by = 1)

Results2$ws10<-seq(from = 0.5, to = 39.5, by = 1)

Results2$ws10_PDF_all<-table(cut(points_processed$wmax_grd,BREAKS, right=FALSE)) / length(points_processed$wmax_grd) 

ws10_PDF_MTS<-points_processed$wmax_grd[-which(points_processed$MTS=='no')]

Results2$ws10_PDF_MTS<-table(cut(ws10_PDF_MTS,BREAKS, right=FALSE)) / length(ws10_PDF_MTS) 


Fig_1d<-ggplot(Results2,aes(x =ws10))+
	geom_path( aes(y = ws10_PDF_all,color = "All"))+
	geom_line( aes(y = ws10_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,40), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 40, by = 5),
	                   labels = every_nth(seq(from = 0, to = 40, by = 5),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.15), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.15, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.15, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Max. Sustained Wind Speed [m/s]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(d)"))+
	theme(axis.text = element_text(colour = "black"))



# for max precipitation   -------------------------------------------------------------

Results3<-data.frame(NO=seq(1,50,1))

BREAKS<-seq(from = 0.00, to = 50, by = 1)

Results3$prmax<-seq(from = 0.5, to = 49.5, by = 1)

Results3$prmax_PDF_all<-table(cut(points_processed$prmax_grd,BREAKS, right=FALSE)) / length(points_processed$prmax_grd) 

prmax_PDF_MTS<-points_processed$prmax_grd[-which(points_processed$MTS=='no')]

Results3$prmax_PDF_MTS<-table(cut(prmax_PDF_MTS,BREAKS, right=FALSE)) / length(prmax_PDF_MTS) 


Fig_1e<-ggplot(Results3,aes(x =prmax))+
	geom_path( aes(y = prmax_PDF_all,color = "All"))+
	geom_line( aes(y = prmax_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,50), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 50, by = 5),
	                   labels = every_nth(seq(from = 0, to = 50, by = 5),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.1), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.1, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.1, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Max. Precipitation [mm]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(e)"))+
	theme(axis.text = element_text(colour = "black"))



# for sum precipitation   -------------------------------------------------------------

Results4<-data.frame(NO=seq(1,32,1))

BREAKS<-seq(from = 0.00, to = 800, by = 25)

Results4$prsum<-seq(from = 12.5, to = 787.5, by = 25)

Results4$prsum_PDF_all<-table(cut(points_processed$prsum_grd,BREAKS, right=FALSE)) / length(points_processed$prsum_grd) 

prsum_PDF_MTS<-points_processed$prsum_grd[-which(points_processed$MTS=='no')]

Results4$prsum_PDF_MTS<-table(cut(prsum_PDF_MTS,BREAKS, right=FALSE)) / length(prsum_PDF_MTS) 


Fig_1f<-ggplot(Results4,aes(x =prsum))+
	geom_path( aes(y = prsum_PDF_all,color = "All"))+
	geom_line( aes(y = prsum_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,800), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 800, by = 50),
	                   labels = every_nth(seq(from = 0, to = 800, by = 50),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.1), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.1, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.1, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('All' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Sum of Precipitation [mm]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(f)"))+
	theme(axis.text = element_text(colour = "black"))


# plotting All  ------------------------------------------
rm_legend <- function(p){p + theme(legend.position = "none")}

commom_legend<- get_legend(Fig_1a)

figure <- ggarrange(rm_legend(Fig_1a),rm_legend(Fig_1b),
	rm_legend(Fig_1c),rm_legend(Fig_1d),
	rm_legend(Fig_1e),rm_legend(Fig_1f),
	ncol = 2, nrow = 3,
	align = c("hv"),
	common.legend = TRUE,
	legend = "bottom",
	legend.grob = commom_legend)


ggsave(plot = figure, filename = paste0("Figure_5-3_.png"),path = ("/home/ahmed/Desktop/FRM_upload/"),
       height = 200, width =175, units = "mm", dpi = 300, device = "png",limitsize = FALSE)


# Initial TC genesis -----------------------------------------------------------

TS_gensis<-points_processed[which(points_processed$Track==1),]


MTS_gensis<-points_processed[which(points_processed$Track==1 & points_processed$MTS=='yes'),]
nonMTS_gensis<-points_processed[which(points_processed$Track==1 & points_processed$MTS=='no'),]

TS_gensis$ID<-seq(1,nrow(TS_gensis),1)

coordinates(TS_gensis) <- ~LON+LAT

proj4string(TS_gensis) <- CRS( "+proj=longlat +datum=WGS84" )

TS_gensis<-st_as_sf(TS_gensis)

# shape file 

coastlines <- st_read('/media/ahmed/Volume/TC_FRM/R/FRM_certficate/continent shapefile/continent.shp')


Fig_2a<- ggplot()+
	geom_sf(data=coastlines, aes(),show.legend = NA,size = 0.1)+

	geom_sf(data=TS_gensis, aes(color=MTS),show.legend = NA,size = 0.1)+
	
	scale_color_manual(values = c('black', "red"))+
	
	scale_x_continuous(limits = c(30,80), expand = c(0, 0)) +
	scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
	
	
	
	theme_bw(base_size = 8)+

	theme(legend.position = 'bottom',
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_blank(),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(5, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(panel.grid.major = element_line(colour = "grey"),
	      panel.background = element_rect(fill = "transparent",colour = NA),
	      panel.ontop=TRUE,
	      plot.background = element_rect(fill = "transparent",colour = NA))+
	
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	labs(title = expression("(a)"))+
	theme(axis.text = element_text(colour = "black"))


# PDFs --------------------------------------------------------------------

Results<-data.frame(NO=seq(1,30,1))

BREAKS<-seq(from = 0.00, to = 30, by = 1)

Results$Lat<-seq(from = 0.5, to = 29.5, by = 1)

Results$Lat_PDF_all<-table(cut(nonMTS_gensis$LAT,BREAKS, right=FALSE)) / length(nonMTS_gensis$LAT) 

Results$Lat_PDF_MTS<-table(cut(MTS_gensis$LAT,BREAKS, right=FALSE)) / length(MTS_gensis$LAT) 

Fig_2b<-ggplot(Results,aes(x =Lat))+
	geom_path( aes(y = Lat_PDF_all,color = "non-MTS"))+
	geom_line( aes(y = Lat_PDF_MTS,color = "MTS"))+
	
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,30), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 30, by = 2),
	                   labels = every_nth(seq(from = 0, to = 30, by = 2),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.15), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.15, by = 0.005),
	                   labels = every_nth(seq(from = 0, to = 0.15, by = 0.005),2, inverse=TRUE))+
	
	scale_colour_manual(values = c('non-MTS' = 'black',
		           'MTS' = 'red'))+
	
	theme(legend.position = 'bottom',
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_blank(),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(5, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(0,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab(expression("Latitude [degree]"))+
	ylab("Relative Frequency")+
	labs(title = expression("(b)"))+
	theme(axis.text = element_text(colour = "black"))

# plotting All  ------------------------------------------
rm_legend <- function(p){p + theme(legend.position = "none")}

commom_legend<- get_legend(Fig_2b)

figure <- ggarrange(rm_legend(Fig_2a),rm_legend(Fig_2b),
	ncol = 2, nrow = 1,
	align = c("hv"),
	common.legend = TRUE,
	legend = "bottom",
	legend.grob = commom_legend)


ggsave(plot = figure, filename = paste0("Figure_5-4_.png"),path = ("/home/ahmed/Desktop/FRM_upload/"),
       height = 100, width =160, units = "mm", dpi = 300, device = "png",limitsize = FALSE)






# END ---------------------------------------------------------------------

