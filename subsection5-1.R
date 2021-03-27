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

# Anomalous---------------------------------------------------------------------
historical<-which(year(ERA5_filtered$time2)<=2005)

Anomalous_df <- list(ERA5_filtered$anomalous,ERA5_filtered$anomalous[historical],CORDEX_filtered$anomalous)

max_length <- max(sapply(Anomalous_df,length))

Anomalous_df<- sapply(Anomalous_df, function(x){
	c(x, rep(NA, max_length - length(x)))
})

Anomalous_df<-data.frame(Anomalous_df)

colnames(Anomalous_df)<-c("ERA5","ERA5-2005","CORDEX")

# Anomalous- plotting ----------------------------------------------------------

# for Relative frequency 

Results<-data.frame(NO=seq(1,80,1))

BREAKS<-seq(from = 0, to = 12, by = 0.15)

Results$Temp<-seq(from = 0.075, to = 11.925, by = 0.15)

Results$ERA5<-table(cut(Anomalous_df$ERA5,BREAKS, right=FALSE)) / length_notNA(Anomalous_df$ERA5) 
Results$`ERA5-2005`<-table(cut(Anomalous_df$`ERA5-2005`,BREAKS, right=FALSE)) / length_notNA(Anomalous_df$`ERA5-2005`) 
Results$CORDEX<-table(cut(Anomalous_df$CORDEX,BREAKS, right=FALSE)) / length_notNA(Anomalous_df$CORDEX) 

Fig_1a<-ggplot(Results,aes(x =Temp))+
	geom_line( aes(y = ERA5,color = "ERA5, 1990-2019"), size = 1)+
	geom_line( aes(y = `ERA5-2005`,color = "ERA5, 1990-2005"))+
	geom_line( aes(y = CORDEX, color = "CORDEX, 1990-2005"))+
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,12), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 12, by = 1),
	                   labels = every_nth(seq(from = 0, to = 12, by = 1),2, inverse=TRUE) ) +
	
	scale_y_continuous(limits = c(0,0.5), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.5, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.5, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('ERA5, 1990-2019' = 'black',
		           'ERA5, 1990-2005' = 'red',
		           "CORDEX, 1990-2005"= "blue"))+
	
	theme(legend.position = "bottom",
	      legend.title = element_blank(),
	      legend.text = element_text(size = 8),
	      legend.box.background = element_blank(),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(5,0,2,0, unit = "mm"))+
	
	theme(plot.margin = margin(5,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	xlab("Vertically Integrated Anomalous Temperature [K]")+
	ylab("Relative Frequency")+
	labs(title = expression("(a)"))+
	theme(axis.text = element_text(colour = "black"))

# plot boxplot 
Anomalous_df_melted<- Anomalous_df

colnames(Anomalous_df_melted)<-c('ERA5, 1990-2019','ERA5, 1990-2005', "CORDEX, 1990-2005")

Anomalous_df_melted<- melt(Anomalous_df_melted)



Fig_1b<-ggplot(data = Anomalous_df_melted, aes(y=value, x = variable, color=variable))+
	
	geom_violin(na.rm = T,scale = "width",lwd=0.3)+
	
	geom_boxplot(width=0.1,outlier.colour = "red",outlier.fill = NA,
	             outlier.size = 0.4,outlier.shape = 1,lwd=0.3)+
	theme_bw(base_size = 8)+
	
	scale_y_continuous(limits = c(0,12), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 12, by = 1),
	                   labels = every_nth(seq(from = 0, to = 12, by = 1),2, inverse=TRUE))+
	
	scale_colour_manual(values = c('black', 'red', "blue"))+
	
	theme(legend.position = c(0.85, 0.85),
	      legend.title = element_blank(),
	      legend.text = element_text(size = 10),
	      legend.box.background = element_rect(colour = "black"),
	      legend.key.width = unit(1, "line"),
	      legend.spacing =  unit(0.05, 'mm'),
	      legend.key.size  = unit(0.15, "cm"),
	      legend.margin = margin(-1,1,0,0, unit = "mm"))+
	
	theme(plot.margin = margin(5,4.5,0,0, unit = "mm"),
	      axis.line = element_line(colour = "black",size=0.25),
	      panel.border = element_rect(colour = "black", fill = NA, size=0.5))+
	
	ylab("Vertically Integrated Anomalous Temperature [K]")+
	labs(title = expression("(b)"))+
	theme(axis.text = element_text(colour = "black"),
	      axis.title.x = element_blank())

# Vorticity---------------------------------------------------------------------
Vorticity_df <- list(ERA5_filtered$vorticity,ERA5_filtered$vorticity[historical],CORDEX_filtered$vorticity)

max_length <- max(sapply(Vorticity_df,length))

Vorticity_df<- sapply(Vorticity_df, function(x){
	c(x, rep(NA, max_length - length(x)))
})

Vorticity_df<-data.frame(Vorticity_df)

colnames(Vorticity_df)<-c("ERA5","ERA5-2005","CORDEX")

# Vorticity- plotting ----------------------------------------------------------

# for Relative frequency 

Results1<-data.frame(NO=seq(1,40,1))

BREAKS<-seq(from = -0.001, to = 0.003, by = 1e-04)

Results1$Vorticity<-seq(from = -0.00095, to = 0.00295, by = 1e-04)

Results1$ERA5<-table(cut(Vorticity_df$ERA5,BREAKS, right=FALSE)) / length_notNA(Vorticity_df$ERA5) 
Results1$`ERA5-2005`<-table(cut(Vorticity_df$`ERA5-2005`,BREAKS, right=FALSE)) / length_notNA(Vorticity_df$`ERA5-2005`) 
Results1$CORDEX<-table(cut(Vorticity_df$CORDEX,BREAKS, right=FALSE)) / length_notNA(Vorticity_df$CORDEX) 

Fig_2a<-ggplot(Results1,aes(x =Vorticity*1e03))+
	geom_path( aes(y = ERA5,color = "ERA5, 1990-2019"))+
	geom_line( aes(y = `ERA5-2005`,color = "ERA5, 1990-2005"))+
	geom_line( aes(y = CORDEX, color = "CORDEX, 1990-2005"))+
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(-2.5,2.5), expand = c(0, 0),
	                   breaks = seq(from = -2.5, to = 2.5, by = 0.25),
	                   labels = every_nth(seq(from = -2.5, to = 2.5, by = 0.25),2, inverse=TRUE)) +
	
	scale_y_continuous(limits = c(0,0.6), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.6, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.6, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('ERA5, 1990-2019' = 'black',
		           'ERA5, 1990-2005' = 'red',
		           "CORDEX, 1990-2005"= "blue"))+
	
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
	
	xlab(expression("Relative Vorticity "(s^{-1}~x10^{3})))+
	ylab("Relative Frequency")+
	labs(title = expression("(c)"))+
	theme(axis.text = element_text(colour = "black"))

# plot boxplot 
Vorticity_df_melted<- Vorticity_df

colnames(Vorticity_df_melted)<-c('ERA5, 1990-2019','ERA5, 1990-2005', "CORDEX, 1990-2005")

Vorticity_df_melted<- melt(Vorticity_df_melted)



Fig_2b<-ggplot(data = Vorticity_df_melted, aes(y=value*1e03, x = variable, color=variable))+
	
	geom_violin(na.rm = T,scale = "width",lwd=0.3)+
	
	geom_boxplot(width=0.1,outlier.colour = "red",outlier.fill = NA,
	             outlier.size = 0.4,outlier.shape = 1,lwd=0.3)+
	theme_bw(base_size = 8)+
	
	scale_y_continuous(limits = c(-1,3), expand = c(0, 0),
	                   breaks = seq(from = -1, to = 3, by = 0.25),
	                   labels = every_nth(seq(from = -1, to = 3, by = 0.25),4, inverse=TRUE))+
	
	scale_colour_manual(values = c('black', 'red', "blue"))+
	
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
	
	ylab(expression("Relative Vorticity "(s^{-1}~x10^{3})))+
	labs(title = expression("(d)"))+
	theme(axis.text = element_text(colour = "black"),
	      axis.title.x = element_blank())

# Wind_speed---------------------------------------------------------------------
Wind_speed_df <- list(ERA5_filtered$surface_wind,ERA5_filtered$surface_wind[historical],CORDEX_filtered$surface_wind)

max_length <- max(sapply(Wind_speed_df,length))

Wind_speed_df<- sapply(Wind_speed_df, function(x){
	c(x, rep(NA, max_length - length(x)))
})

Wind_speed_df<-data.frame(Wind_speed_df)

colnames(Wind_speed_df)<-c("ERA5","ERA5-2005","CORDEX")

# Wind_speed plotting ----------------------------------------------------------

# for Relative frequency 

Results2<-data.frame(NO=seq(1,25,1))

BREAKS<-seq(from = 0, to = 25, by = 1)

Results2$wind_speed<-seq(from = 0.5, to = 24.5, by = 1)

Results2$ERA5<-table(cut(Wind_speed_df$ERA5,BREAKS, right=FALSE)) / length_notNA(Wind_speed_df$ERA5) 
Results2$`ERA5-2005`<-table(cut(Wind_speed_df$`ERA5-2005`,BREAKS, right=FALSE)) / length_notNA(Wind_speed_df$`ERA5-2005`) 
Results2$CORDEX<-table(cut(Wind_speed_df$CORDEX,BREAKS, right=FALSE)) / length_notNA(Wind_speed_df$CORDEX) 

Fig_3a<-ggplot(Results2,aes(x =wind_speed))+
	geom_line( aes(y = ERA5,color = "ERA5, 1990-2019"), size = 1)+
	geom_line( aes(y = `ERA5-2005`,color = "ERA5, 1990-2005"))+
	geom_line( aes(y = CORDEX, color = "CORDEX, 1990-2005"))+
	theme_bw(base_size = 8)+
	scale_x_continuous(limits = c(0,20), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 20, by = 1),
	                   labels = every_nth(seq(from = 0, to = 20, by = 1),2, inverse=TRUE) ) +
	
	scale_y_continuous(limits = c(0,0.25), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 0.25, by = 0.01),
	                   labels = every_nth(seq(from = 0, to = 0.25, by = 0.01),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('ERA5, 1990-2019' = 'black',
		           'ERA5, 1990-2005' = 'red',
		           "CORDEX, 1990-2005"= "blue"))+
	
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
	
	xlab("Wind Speed [m/s]")+
	ylab("Relative Frequency")+
	labs(title = expression("(e)"))+
	theme(axis.text = element_text(colour = "black"))

# plot boxplot 
Wind_speed_df_melted<- Wind_speed_df

colnames(Wind_speed_df_melted)<-c('ERA5, 1990-2019','ERA5, 1990-2005', "CORDEX, 1990-2005")

Wind_speed_df_melted<- melt(Wind_speed_df_melted)



Fig_3b<-ggplot(data = Wind_speed_df_melted, aes(y=value, x = variable, color=variable))+
	
	geom_violin(na.rm = T,scale = "width",lwd=0.3)+
	
	geom_boxplot(width=0.1,outlier.colour = "red",outlier.fill = NA,
	             outlier.size = 0.4,outlier.shape = 1,lwd=0.3)+
	theme_bw(base_size = 8)+
	
	scale_y_continuous(limits = c(0,25), expand = c(0, 0),
	                   breaks = seq(from = 0, to = 25, by = 1),
	                   labels = every_nth(seq(from = 0, to = 25, by = 1),5, inverse=TRUE))+
	
	scale_colour_manual(values = c('black', 'red', "blue"))+
	
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
	
	ylab("Wind Speed [m/s]")+
	labs(title = expression("(f)"))+
	theme(axis.text = element_text(colour = "black"),
	      axis.title.x = element_blank())









# plotting All  ------------------------------------------
rm_legend <- function(p){p + theme(legend.position = "none")}

commom_legend<- get_legend(Fig_1a)

figure <- ggarrange(rm_legend(Fig_1a),rm_legend(Fig_1b),
	rm_legend(Fig_2a),rm_legend(Fig_2b),
	rm_legend(Fig_3a),rm_legend(Fig_3b),
	ncol = 2, nrow = 3,
	align = c("hv"),
	common.legend = TRUE,
	legend = "bottom",
	legend.grob = commom_legend)


ggsave(plot = figure, filename = paste0("Figure_5-1_.png"),path = ("/home/ahmed/Desktop/FRM_upload/"),
       height = 200, width =175, units = "mm", dpi = 300, device = "png",limitsize = FALSE)


# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
	if( equal.variance==FALSE ) 
	{
		se <- sqrt( (s1^2/n1) + (s2^2/n2) )
		# welch-satterthwaite df
		df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
	} else
	{
		# pooled standard deviation, scaled by the sample sizes
		se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
		df <- n1+n2-2
	}      
	t <- (m1-m2-m0)/se 
	dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
	names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
	return(dat) 
}

m1<-mean(Anomalous_df$ERA5,na.rm = T)
m2<-mean(Anomalous_df$`ERA5-2005`,na.rm = T)

s1<-sd(Anomalous_df$ERA5,na.rm = T)
s2<-sd(Anomalous_df$`ERA5-2005`,na.rm = T)

n1<-length_notNA(Anomalous_df$ERA5)
n2<-length_notNA(Anomalous_df$`ERA5-2005`)


tt2<-t.test2(m1,m2,s1,s2,n1,n2)

tt <-t.test(Anomalous_df$ERA5,Anomalous_df$`ERA5-2005`,var.equal = T)

tt$statistic == tt2[["t"]]

tt$p.value == tt2[["p-value"]]


var.test(Anomalous_df$ERA5,Anomalous_df$`ERA5-2005`,alternative = "two.sided")

t.test(ERA5_filtered$anomalous[-historical],ERA5_filtered$anomalous[historical])
