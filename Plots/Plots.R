
setwd("C:Folder/Plots") # Choose Directory

library("ggplot2")
library(lubridate)
library("ggthemes")
library(latex2exp)
library(egg) 


ramp <- colorRamp(c("green", "blue", "navy"))
ramp.list <- rgb( ramp(seq(0, 1, length = 12)), max = 255)


## Getting data
plot<-read.csv("RISK.TAIL.csv", sep = ";")
head(plot)
plot<-plot[,2:5]
head(plot)

## Gettin date 
datum<-read.csv("datum.csv", sep=";")
head(datum)

dim(plot)[2]
dim(datum)[1]
tail(plot,100)
fechas <- as.POSIXct(datum[(1:dim(datum)[1]),1], tz="",format = "%d-%m-%Y")
fechas


D1<-data.frame(fechas, plot[1:dim(datum)[1]                       ,1:dim(plot)[2]])
D2<-data.frame(fechas, plot[(1*dim(datum)[1]+1):(2*dim(datum)[1]) ,1:dim(plot)[2]])
D3<-data.frame(fechas, plot[(2*dim(datum)[1]+1):(3*dim(datum)[1]) ,1:dim(plot)[2]])
D4<-data.frame(fechas, plot[(3*dim(datum)[1]+1):(4*dim(datum)[1]) ,1:dim(plot)[2]])
D5<-data.frame(fechas, plot[(4*dim(datum)[1]+1):(5*dim(datum)[1]) ,1:dim(plot)[2]])


## Plotting individually

D1<-ggplot(D1[50:dim(D1)[1],], aes(x=fechas, y=retornos))+
  geom_point(aes(x=fechas, y = retornos),col="gray20", size=0.2)+
  geom_line(aes(x=fechas, y =  var), col= "green")+
  geom_line(aes(x=fechas, y = es), col="blue")+
  theme(panel.spacing = unit(0.5,"lines"))+
  theme(panel.background = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.line=element_line(color="gray80"))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.position="none")+
  labs(x = "",  y = "APPL")+
  scale_y_continuous(limits = c(0.01, 0.12), breaks=c(0.04,0.08,0.12))+
  scale_x_datetime(breaks = ("2 year"), date_labels = c("%Y"))+ 
  theme(axis.text=element_text(size=12))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  
  theme(axis.ticks.y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
D1

D2<-ggplot(D2[50:dim(D2)[1],], aes(x=fechas, y=retornos))+
  geom_point(aes(x=fechas, y = retornos),col="gray20", size=0.2)+
  geom_line(aes(x=fechas, y =  var), col= "green")+
  geom_line(aes(x=fechas, y = es), col="blue")+
  theme(panel.spacing = unit(0.5,"lines"))+
  theme(panel.background = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.line=element_line(color="gray80"))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.position="none")+
  labs(x = "",  y = "AMZN")+
  scale_y_continuous(limits = c(0.01, 0.12), breaks=c(0.04,0.08,0.12))+
  scale_x_datetime(breaks = ("2 year"), date_labels = c("%Y"))+ 
  theme(axis.text=element_text(size=12))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
D2

D3<-ggplot(D3[50:dim(D3)[1],], aes(x=fechas, y=retornos))+
  geom_point(aes(x=fechas, y = retornos),col="gray20", size=0.2)+
  geom_line(aes(x=fechas, y =  var), col= "green")+
  geom_line(aes(x=fechas, y = es), col="blue")+
  theme(panel.spacing = unit(0.5,"lines"))+
  theme(panel.background = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.line=element_line(color="gray80"))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.position="none")+
  labs(x = "",  y = "CSCO")+
  scale_y_continuous(limits = c(0.01, 0.12), breaks=c(0.04,0.08,0.12))+
  scale_x_datetime(breaks = ("2 year"), date_labels = c("%Y"))+ 
  theme(axis.text=element_text(size=12))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
D3

D4<-ggplot(D4[50:dim(D4)[1],], aes(x=fechas, y=retornos))+
  geom_point(aes(x=fechas, y = retornos),col="gray20", size=0.2)+
  geom_line(aes(x=fechas, y =  var), col= "green")+
  geom_line(aes(x=fechas, y = es), col="blue")+
  theme(panel.spacing = unit(0.5,"lines"))+
  theme(panel.background = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.line=element_line(color="gray80"))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.position="none")+
  labs(x = "",  y = "IBM")+
  scale_y_continuous(limits = c(0.01, 0.12), breaks=c(0.04,0.08,0.12))+
  scale_x_datetime(breaks = ("2 year"), date_labels = c("%Y"))+ 
  theme(axis.text=element_text(size=12))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
D4

D5<-ggplot(D5[50:dim(D5)[1],], aes(x=fechas, y=retornos))+
  geom_point(aes(x=fechas, y = retornos),col="gray20", size=0.3)+
  geom_line(aes(x=fechas, y =  var), col= "green")+
  geom_line(aes(x=fechas, y = es), col="blue")+
  theme(panel.spacing = unit(0.5,"lines"))+
  theme(panel.background = element_blank())+
  theme(strip.background = element_blank())+
  theme(axis.line=element_line(color="gray80"))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.position="none")+
  labs(x = "",  y = "MSFT")+
  scale_y_continuous(limits = c(0.01, 0.12), breaks=c(0.04,0.08,0.12))+
  scale_x_datetime(breaks = ("1 year"), date_labels = c("%Y"))+ 
  theme(axis.text=element_text(size=12))+
  theme(axis.ticks.y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
D5

## Plotting all together
P<-ggarrange(D1,D2,D3,D4,D5,ncol = 1)
P
















