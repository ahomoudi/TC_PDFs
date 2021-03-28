library(ff)
library(ffbase)
library(data.table)
library(stringr)
library(reshape2)
library(viridis)
library(tidyverse)

library(ggplot2)
library(ggpubr)
library(scales)

library(reshape2)
library(bivariate)
library(dplyr)
library(fields)

source("/media/ahmed/Daten/GCMs/ITCZ_analysis/extra_function.R")
# read ERA5 warm cores points -------------------------------------------------

ERA5_filtered <- read.csv("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/all_filtered.csv")

ERA5_filtered <- ERA5_filtered[,-c(1,2)]


# check that all points in the tropical cyclone season
library(PCICt)
library(lubridate)

time.origin <- PCICt::as.PCICt.default("1900-01-01 00:00:00.0", 
	                   cal = "gregorian")

ERA5_filtered$time2<- PCICt::as.PCICt(ERA5_filtered $TIME*3600, 
	                  origin = time.origin,cal = "gregorian")

ERA5_filtered$month<- month(PCICt::as.PCICt(ERA5_filtered $TIME*3600, 
		    origin = time.origin,cal = "gregorian"))

plot(ERA5_filtered$month)

# read CORDEX warm cores points ------------------------------------------------
CORDEX_filtered <- read.csv("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/all_filtered_CORDEX.csv")

CORDEX_filtered <- CORDEX_filtered[,-c(1)]

# check that all points in the tropical cyclone season

time.origin <- PCICt::as.PCICt.default("1949-12-01 00:00:00", 
	                   cal = "365_day")

CORDEX_filtered$time2<- PCICt::as.PCICt(CORDEX_filtered$TIME*3600*24, 
		origin = time.origin,cal = "365_day")

CORDEX_filtered$month<- month(PCICt::as.PCICt(CORDEX_filtered$TIME*3600*24, 
		      origin = time.origin,cal = "365_day"))

nop<-which(CORDEX_filtered$month!=c(5,6,9,10,11))

plot(CORDEX_filtered$month)

CORDEX_filtered<-CORDEX_filtered[-nop,]

plot(CORDEX_filtered$month)

length_notNA<-function(x){
	length(x[!is.na(x)])
}


# joint PDF  ERA5 #######################################

# ERA5 1  -----------------------------------------------------------------


df_ERA5_1<- ERA5_filtered %>% dplyr::select("vorticity","anomalous")

colnames(df_ERA5_1)<-c("x","y")

df_ERA5_1$x<-df_ERA5_1$x/sd(df_ERA5_1$x)

df_ERA5_1$y<-df_ERA5_1$y/sd(df_ERA5_1$y)


dat1 <- as.data.frame(df_ERA5_1)
f1 = kbvpdf (dat1$x, dat1$y, 5, 5)
plot (f1, TRUE, xlab="vorticity", ylab="anomalous")
x1<- bvmat(f1,n=50)

PDFs1<-matrix(x1[["fv"]],nrow = 50, byrow = F)
PDFs1<-array(PDFs1, dim = c(50,50), dimnames = list(vorticity = x1[["x"]] , anomalous = x1[["y"]] ))

# ERA5 2  -----------------------------------------------------------------

df_ERA5_2<- ERA5_filtered %>% dplyr::select("vorticity","surface_wind")

colnames(df_ERA5_2)<-c("x","y")

df_ERA5_2$x<-df_ERA5_2$x/sd(df_ERA5_2$x)

df_ERA5_2$y<-df_ERA5_2$y/sd(df_ERA5_2$y)


dat2 <- as.data.frame(df_ERA5_2)
f2 = kbvpdf (dat2$x, dat2$y, 5, 5)
plot (f2, TRUE, xlab="vorticity", ylab="surface_wind")
x2<- bvmat(f2,n=50)

PDFs2<-matrix(x2[["fv"]],nrow = 50, byrow = F)
PDFs2<-array(PDFs2, dim = c(50,50), dimnames = list(vorticity = x2[["x"]] , surface_wind = x2[["y"]] ))

#Divide the screen in 1 columns and 2 lines
par(mfrow=c(1,2))


image.plot(PDFs1)

image.plot(PDFs2)

# joint PDF  CORDEX #######################################

# CORDEX 1  -----------------------------------------------------------------


df_CORDEX_1<- CORDEX_filtered %>% dplyr::select("vorticity","anomalous")

colnames(df_CORDEX_1)<-c("x","y")

df_CORDEX_1$x<-df_CORDEX_1$x/sd(df_CORDEX_1$x)

df_CORDEX_1$y<-df_CORDEX_1$y/sd(df_CORDEX_1$y)


dat3 <- as.data.frame(df_CORDEX_1)
f3= kbvpdf (dat3$x, dat3$y, 5, 5)
plot (f3, TRUE, xlab="vorticity", ylab="anomalous")
x3<- bvmat(f3,n=50)

PDFs3<-matrix(x3[["fv"]],nrow = 50, byrow = F)
PDFs3<-array(PDFs3, dim = c(50,50), dimnames = list(vorticity = x3[["x"]] , anomalous = x3[["y"]] ))

# CORDEX 2  -----------------------------------------------------------------

df_CORDEX_2<- CORDEX_filtered %>% dplyr::select("vorticity","surface_wind")

colnames(df_CORDEX_2)<-c("x","y")

df_CORDEX_2$x<-df_CORDEX_2$x/sd(df_CORDEX_2$x)

df_CORDEX_2$y<-df_CORDEX_2$y/sd(df_CORDEX_2$y)


dat4 <- as.data.frame(df_CORDEX_2)
f4 = kbvpdf (dat4$x, dat4$y, 5, 5)
plot (f4, TRUE, xlab="vorticity", ylab="surface_wind")
x4<- bvmat(f4,n=50)

PDFs4<-matrix(x4[["fv"]],nrow = 50, byrow = F)
PDFs4<-array(PDFs4, dim = c(50,50), dimnames = list(vorticity = x4[["x"]] , surface_wind = x4[["y"]] ))

#Divide the screen in 2 columns and 2 lines
par(mfrow=c(2,2))


image.plot(PDFs1)

image.plot(PDFs2)

image.plot(PDFs3)

image.plot(PDFs4)

# plotting ----------------------------------------------------------------
library(reshape2)
library(scales)
library(viridis)

jet.colors <- colorRampPalette(c("white", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


# Fig a -------------------------------------------------------------------

Dfa<-reshape2::melt(PDFs1)

fig_a<- ggplot()+ 
	geom_raster(data = Dfa , aes(x = vorticity, y = anomalous, fill = value))+ 
	
	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), breaks = seq(0,0.006,0.0005),oob=squish,
		 na.value="white",
		 guide = guide_colourbar(direction = "horizontal",
		 	    title.position = "left",
		 	    label.position = "bottom"))+
	scale_y_continuous(expand = c(0,0), limits = c(-5,40), breaks = seq(-5,40,5))+
	scale_x_continuous(expand = c(0,0), limits = c(-15,35), breaks = seq(-15,35,5))+
	theme_light(base_size = 11)+
	
	theme(panel.grid = element_blank(),
	      panel.background = element_rect(fill = "transparent",colour = NA),
	      panel.ontop=TRUE,
	      plot.background = element_rect(fill = "transparent",colour = NA),
	      plot.margin = unit(c(5,5,5,5), "mm"),
	      plot.title = element_text(hjust = 0,margin =margin (0,0,5,0),lineheight = 0.1),
	      legend.position = "bottom",
	      legend.key.width =unit(30.0, 'mm'),
	      legend.key.height =unit(5.0, 'mm'),
	      legend.title = element_text(hjust = 0.65))+
	guides(fill = guide_colourbar(title.position = "bottom",title.hjust = .3,label.position = "bottom"))+
	theme(axis.title = element_text(colour = "black", size = 10),
	      axis.text = element_text(colour = "black", size = 6))+
	labs(fill= expression(~~~~~" Joint Relative Frequency"),
	     title = expression("(a)"))+
	ylab(expression(paste("T/",sigma[T])))+
	xlab(expression(paste(zeta,"/",sigma[zeta])))

# Fig b -------------------------------------------------------------------
Dfb<-reshape2::melt(PDFs2)

fig_b<- ggplot() +
	geom_raster(data = Dfb , aes(x = vorticity, y = surface_wind, fill = value)) + 
	
	coord_quickmap()+
	
	
	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), 
		 breaks = seq(0,0.006,0.0005),oob=squish,
		 na.value="white",
		 guide = guide_colourbar(direction = "horizontal",
		 	    title.position = "left",
		 	    label.position = "bottom"))+
	
	scale_y_continuous(expand = c(0,0), limits = c(-5,15), breaks = seq(-5,15,5))+
	scale_x_continuous(expand = c(0,0), limits = c(-15,35), breaks = seq(-15,35,5))+
	theme_light(base_size = 11) +
	
	theme(panel.grid = element_blank(),
	      panel.background = element_rect(fill = "transparent",colour = NA),
	      panel.ontop=TRUE,
	      plot.background = element_rect(fill = "transparent",colour = NA),
	      plot.margin = unit(c(5,5,5,5), "mm"),
	      plot.title = element_text(hjust = 0,margin =margin (0,0,5,0),lineheight = 0.1),
	      legend.position = "bottom",
	      legend.key.width =unit(25.0, 'mm'),
	      legend.key.height =unit(5.0, 'mm'),
	      legend.title = element_text(hjust = 0.65))+
	guides(fill = guide_colourbar(title.position = "bottom",title.hjust = .3,label.position = "bottom"))+
	theme(axis.title = element_text(colour = "black", size = 10),
	      axis.text = element_text(colour = "black", size = 6))+
	labs(fill= "    Joint Relative Frequency",
	     title = expression("(b)"))+
	ylab(expression(paste("v/",sigma[v])))+
	xlab(expression(paste(zeta,"/",sigma[zeta])))

# Fig c -------------------------------------------------------------------
Dfc<-reshape2::melt(PDFs3)

fig_c<- ggplot()+ 
	geom_raster(data = Dfc , aes(x = vorticity, y = anomalous, fill = value))+ 
	
	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), breaks = seq(0,0.006,0.0005),oob=squish,
		 na.value="white",
		 guide = guide_colourbar(direction = "horizontal",
		 	    title.position = "left",
		 	    label.position = "bottom"))+
	scale_y_continuous(expand = c(0,0), limits = c(-5,40), breaks = seq(-5,40,5))+
	scale_x_continuous(expand = c(0,0), limits = c(-15,35), breaks = seq(-15,35,5))+
	theme_light(base_size = 11)+
	
	theme(panel.grid = element_blank(),
	      panel.background = element_rect(fill = "transparent",colour = NA),
	      panel.ontop=TRUE,
	      plot.background = element_rect(fill = "transparent",colour = NA),
	      plot.margin = unit(c(5,5,5,5), "mm"),
	      plot.title = element_text(hjust = 0,margin =margin (0,0,5,0),lineheight = 0.1),
	      legend.position = "bottom",
	      legend.key.width =unit(25.0, 'mm'),
	      legend.key.height =unit(5.0, 'mm'),
	      legend.title = element_text(hjust = 0.65))+
	guides(fill = guide_colourbar(title.position = "bottom",title.hjust = .3,label.position = "bottom"))+
	theme(axis.title = element_text(colour = "black", size = 10),
	      axis.text = element_text(colour = "black", size = 6))+
	labs(fill= expression(~~~~~" Joint Relative Frequency"),
	     title = expression("(c)"))+
	ylab(expression(paste("T/",sigma[T])))+
	xlab(expression(paste(zeta,"/",sigma[zeta])))
# Fig d -------------------------------------------------------------------
Dfd<-reshape2::melt(PDFs4)

fig_d<- ggplot() +
	geom_raster(data = Dfd , aes(x = vorticity, y = surface_wind, fill = value)) + 
	
	coord_quickmap()+
	
	
	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), 
		 breaks = seq(0,0.006,0.0005),oob=squish,
		 na.value="white",
		 guide = guide_colourbar(direction = "horizontal",
		 	    title.position = "left",
		 	    label.position = "bottom"))+
	
	scale_y_continuous(expand = c(0,0), limits = c(-5,15), breaks = seq(-5,15,5))+
	scale_x_continuous(expand = c(0,0), limits = c(-15,35), breaks = seq(-15,35,5))+
	theme_light(base_size = 11) +
	
	theme(panel.grid = element_blank(),
	      panel.background = element_rect(fill = "transparent",colour = NA),
	      panel.ontop=TRUE,
	      plot.background = element_rect(fill = "transparent",colour = NA),
	      plot.margin = unit(c(5,5,5,5), "mm"),
	      plot.title = element_text(hjust = 0,margin =margin (0,0,5,0),lineheight = 0.1),
	      legend.position = "bottom",
	      legend.key.width =unit(25.0, 'mm'),
	      legend.key.height =unit(5.0, 'mm'),
	      legend.title = element_text(hjust = 0.65))+
	guides(fill = guide_colourbar(title.position = "bottom",title.hjust = .3,label.position = "bottom"))+
	theme(axis.title = element_text(colour = "black", size = 10),
	      axis.text = element_text(colour = "black", size = 6))+
	labs(fill= "    Joint Relative Frequency",
	     title = expression("(d)"))+
	ylab(expression(paste("v/",sigma[v])))+
	xlab(expression(paste(zeta,"/",sigma[zeta])))
# All -------------------------------------------------------------------
library(ggpubr)
rm_legend <- function(p){p + theme(legend.position = "none")}

commom_legend<- get_legend(fig_a)

figure <- ggarrange(rm_legend(fig_a),rm_legend(fig_b),
	rm_legend(fig_c),rm_legend(fig_d),
	ncol = 2, nrow = 2,
	align = c("hv"),
	common.legend = TRUE,
	legend = "bottom",
	legend.grob = commom_legend)


ggsave(plot = figure, filename = paste0("Figure_5-2_.png"),path = ("/home/ahmed/Desktop/FRM_upload/"),
       height = 200, width =175, units = "mm", dpi = 300, device = "png",limitsize = FALSE)




# End ---------------------------------------------------------------------


