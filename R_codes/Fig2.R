rm(list=ls())
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
plot1 <- readPNG('../Plots/Fig2AB.png')
plot2 <- readPNG('../Plots/Fig2CD.png')
plot3 <- readPNG('../Plots/Fig2EF.png')
plot4 <- readPNG('../Plots/Fig2GH.png')


tmp <- arrangeGrob(rasterGrob(plot1),arrangeGrob(rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),ncol=3),ncol=1,heights=c(10,8))
ggsave('../Plots/Fig2.png',tmp,width=8,height=10)