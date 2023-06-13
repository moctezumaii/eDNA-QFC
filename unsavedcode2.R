#Spatial corrrelation using new distances
#Get distances from the "working code"
library(gtools)
library(doBy)
library(reshape)
library(ggplot2)
library(gdistance)
library(tidyverse)
library(raster)
library(fasterize)
library(terra)
library(vegan)
library(nlme)
library(raster)
library(fasterize)
library(sf)
library(rgdal)


### NEW CODE TO GET CORRELATIONS

for(i in unique(subDat$lake)){
  
  #  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
  par(mfrow=c(4,2))
  z<-0
  for(t in 1:4){  
    print(paste0(i," ",t))
    #allDat2<-subDat%>%filter(gear == t)
    sr_dist<-readRDS(paste0(substr(i,1,9),"dist_xy_4.RDS"))
    
    Datanew<-allDat2%>%filter(lake == i)
    Datasp<-Datanew[,8:ncol(Datanew)]
    if(nrow(Datasp)<2){
      print(paste0(i, " did not have enough data for", t))
      z<-z+1
      next
      
    }
    #print(paste0(i, " had a z of ", z))
    rowstokeep<-which(rowSums(Datasp)!=0)
    data.xy<-Datanew[,4:5]
    #data.xy<-data.xy[rowstokeep,]
    break
    #Datasp<-Datasp%>%select(which(!colSums(Datasp, na.rm=TRUE) %in% 0))
    Datasp<-Datasp[rowstokeep,]
    sr_dist<-dist(Datasp)
    xy_dist<-dist(data.xy)
    #  plot(xy_dist, sr_dist)
    #  abline(lm(sr_dist ~ xy_dist), lwd=3, col='red')
    #  lines(lowess(xy_dist, sr_dist), lwd=3, col='pink')
    comm_dist <- vegdist(Datasp)
    #  plot(xy_dist, comm_dist)
    #  abline(lm(comm_dist ~ xy_dist), lwd=3, col='red')
    #  lines(lowess(xy_dist, comm_dist), lwd=3, col='pink')
    #  lines(lowess(xy_dist, comm_dist, f=0.1), lwd=3, col='blue')
    comm_mantel <- mantel(xy_dist, comm_dist)
    #  print(comm_mantel)
    sr_corlog <- mantel.correlog(sr_dist, xy_dist)
    comm_corlog <- mantel.correlog(comm_dist, xy_dist)
    #  print(sr_corlog)
    #  print(comm_corlog)
    #par(mfrow=c(4,2))
    plot(sr_corlog)
    mtext(side=3, paste0('Species Richness ', t))
    plot(comm_corlog)
    mtext(side=3, paste0('Community Composition ', t))
    
    
    break
  }
  mtext(i, side = 3, line = -2, outer = TRUE)
  break
  #  dev.off()
}