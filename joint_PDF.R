library(ggplot2)


#all_filtered <- read.csv("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/all_filtered.csv")
all_filtered <- read.csv("/media/ahmed/Volume/TC_FRM/R/TC_PDFs/all_filtered_CORDEX.csv")

p<-ggplot(all_filtered, aes(x=vorticity)) + 
	geom_histogram(color="black", fill="white", binwidth = 0.00001)
p

length_notNA<-function(x){
	length(x[!is.na(x)])
}

# anomalous#################################################################
x<-sd(all_filtered$anomalous)


max(all_filtered$anomalous/x);min(all_filtered$anomalous/x)


Result_PDF<-data.frame(NO=seq(1,99,1))

BREAKS<-seq(from = 40, to = 0, by = -0.401)

Result_PDF$temp<-seq(from = 39.599, to = 0, by = -0.401)

Result_PDF$temp_PDF<-table(cut(c(all_filtered$anomalous/x),BREAKS, right=FALSE)) / length_notNA(c(all_filtered$anomalous/x)) 

# vorticity #######################################
y<-sd(all_filtered$vorticity)

max(all_filtered$vorticity/y);min(all_filtered$vorticity/y)


BREAKS<-seq(from = 30, to = -10, by = -0.401)

Result_PDF$vorticity<-seq(from = 29.8, to = -9.8, by = -0.401)

Result_PDF$vorticity_PDF<-table(cut(c(all_filtered$vorticity/y),BREAKS, right=FALSE)) / length_notNA(c(all_filtered$vorticity/y)) 

Result_PDF$joint_PDF<-Result_PDF$temp_PDF * Result_PDF$vorticity_PDF

# joint PDF 1 #######################################
library(dplyr)
df<- all_filtered %>% dplyr::select("vorticity","anomalous")

colnames(df)<-c("x","y")

df$x<-df$x/sd(df$x)

df$y<-df$y/sd(df$y)

library(bivariate)


dat <- as.data.frame(df)
f = kbvpdf (dat$x, dat$y, 5, 5)
plot (f, TRUE, xlab="vorticity", ylab="anomalous")
x<- bvmat(f,n=50)

PDFs<-matrix(x[["fv"]],nrow = 50, byrow = F)
PDFs<-array(PDFs, dim = c(50,50), dimnames = list(vorticity = x[["x"]] , anomalous = x[["y"]] ))
library(fields)

image.plot(PDFs)
# joint PDF ###################################################
df<- all_filtered %>% dplyr::select("vorticity","surface_wind")

colnames(df)<-c("x","y")

df$x<-df$x/sd(df$x)

df$y<-df$y/sd(df$y)

library(bivariate)


dat <- as.data.frame(df)
f1 = kbvpdf (dat$x, dat$y,5,5)
plot (f1, TRUE, xlab="vorticity", ylab="surface_wind")

x1<- bvmat(f1,n=50)

PDFs1<-matrix(x1[["fv"]],nrow = 50, byrow = F)
PDFs1<-array(PDFs1, dim = c(50,50), dimnames = list(vorticity = x1[["x"]] , surface_wind = x1[["y"]] ))

image.plot(PDFs1)

#plotting#####################################
library(dplyr)
library(contoureR)  
library(raster)
library(reshape2)
library(scales)
library(viridis)

# fig a ------------------------------------------------------------------
Df<-reshape2::melt(PDFs)

jet.colors <- colorRampPalette(c("white", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

result<-as.data.frame(getContourLines(x=Df$vorticity,
	                  y=Df$anomalous,
	                  z=Df$value,
	                  binwidth = 0.002), stringsAsFactors=TRUE)
             
fig_a<- ggplot() +
	geom_raster(data = Df , aes(x = vorticity, y = anomalous, fill = value)) + 
	
	
	#geom_path(data=result,aes(x,y,group=Group,colour=z)) +
	
             	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), breaks = seq(0,0.006,0.0005),oob=squish,
             		 na.value="white",
             		 guide = guide_colourbar(direction = "horizontal",
             		 	    title.position = "left",
             		 	    label.position = "bottom"))+
	scale_y_continuous(expand = c(0,0),limits = c(-5,15), breaks = seq(-5,15,1))+
	#scale_x_continuous(expand = c(0,0),limits = c(0,15), breaks = seq(0,15,1))+
	scale_x_continuous(expand = c(0,0),limits = c(0,10), breaks = seq(0,10,1))+
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
	labs(fill= expression(~~~~~" Joint Relative Frequency"),
	     title = expression("(a)"))+
	ylab(expression(paste("T/",sigma[T])))+
	xlab(expression(paste(zeta,"/",sigma[zeta])))

# fig b  ------------------------------------------------------------------


Df1<-reshape2::melt(PDFs1)

#jet.colors <- colorRampPalette(c("white", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

result1<-as.data.frame(getContourLines(x=Df1$vorticity,
	                  y=Df1$surface_wind,
	                  z=Df1$value,
	                  binwidth = 0.002), stringsAsFactors=TRUE)

fig_b<- ggplot() +
	geom_raster(data = Df1 , aes(x = vorticity, y = surface_wind, fill = value)) + 
	
	coord_quickmap()+
	
	#geom_path(data=result,aes(x,y,group=Group,colour=z)) +
	
	scale_fill_gradientn(colors = jet.colors(7),limits = c(0,0.006), 
		 breaks = seq(0,0.006,0.0005),oob=squish,
		 na.value="white",
		 guide = guide_colourbar(direction = "horizontal",
		 	    title.position = "left",
		 	    label.position = "bottom"))+
	#scale_y_continuous(expand = c(0,0),limits = c(-5,15), breaks = seq(-5,15,1))+
	#scale_x_continuous(expand = c(0,0),limits = c(0,15), breaks = seq(0,15,1))+
	
	scale_y_continuous(expand = c(0,0),limits = c(-5,12), breaks = seq(-5,12,1))+
	scale_x_continuous(expand = c(0,0),limits = c(0,10), breaks = seq(0,10,1))+
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


# All ---------------------------------------------------------------------
library(ggpubr)
rm_legend <- function(p){p + theme(legend.position = "none")}

commom_legend<- get_legend(fig_a)

ffigure_ALL<- ggarrange(rm_legend(fig_a),rm_legend(fig_b),
	        ncol = 2,nrow = 1,
	        legend = "bottom", align = "hv",
	        common.legend = TRUE,
	        legend.grob = commom_legend)

#ggsave(plot = ffigure_ALL, filename = paste0("Joint_PDF.png"),
       #height = 100, width =175, units = "mm", dpi = 300, device = "png",limitsize = FALSE)             	

ggsave(plot = ffigure_ALL, filename = paste0("Joint_PDF_CORDEX.png"),
       height = 100, width =175, units = "mm", dpi = 300, device = "png",limitsize = FALSE)   
