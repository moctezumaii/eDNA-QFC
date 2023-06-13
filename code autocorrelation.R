
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

#################################
#####     Read in Data        ###
#################################
dat12s<-read.csv("data/12S_100m_22lakes_eDNA_trad_combined_only_spp_level_110719.csv")
dat16s<-read.csv("data/16S_100m_22lakes_eDNA_trad_combined_only_spp_level_110719.csv")
head(dat12s)
#add gene to gear
dat12s$gear<-factor(dat12s$gear)
dat12s$gear<-factor(dat12s$gear,c(levels(dat12s$gear),"benthic_12s","surface_12s"))
dat12s$gear[which(dat12s$gear%in%c("benthic","surface"))]<-paste0(dat12s$gear[which(dat12s$gear%in%c("benthic","surface"))],"_12s")
dat16s$gear<-factor(dat16s$gear)
dat16s$gear<-factor(dat16s$gear,c(levels(dat16s$gear),"benthic_16s","surface_16s"))
dat16s$gear[which(dat16s$gear%in%c("benthic","surface"))]<-paste0(dat16s$gear[which(dat16s$gear%in%c("benthic","surface"))],"_16s")
allDat<-droplevels(smartbind(dat12s,dat16s))
rownames(allDat) <- c()
#reorgnize
allDat<-allDat[,c(1:5,80,81,6:79,82:length(allDat))]
#convert NAs to 0
allDat[,8:length(allDat)][is.na(allDat[,8:length(allDat)])] <- 0
#remove duplicates-- traditional data is duplicated
subDat<-unique(allDat)
#fix brevort lake typo (convert to Brevoort)
subDat$lake[which(subDat$lake=="Brevort Lake")]<-"Brevoort Lake"
subDat<-droplevels(subDat)
#fix gear labels-- change to all lowercase
subDat$gear<-as.factor(tolower(subDat$gear))
allDat
library(vegan)
library(nlme)
####TEST
data(mite)
mite
data(mite.env)
data(mite.xy)
mite.env
sr <- rowSums(mite > 0)
sr
sr_dist <- dist(sr)
sr_dist
data(mite.xy)
mite.xy
########
alldat
allDat
########


data.xy<-subDat[,4:5]
data.xy
allDat
#sr_dist<-dist(data)
#data
########
Datanew<-subDat[,8:ncol(subDat)]
sr_dist<-dist(Datanew)
xy_dist<-dist(data.xy)
xy_dist
#plot(xy_dist, sr_dist)

########
unique(allDat$gear)
allDat2<-subDat%>%filter(gear == "surface_12s")
allDat2
mite.xy
#allDat2<-allDat%>%filter(gear == "surface_12s")

unique(subDat$gear)

geartotest<-c("surface_12s", "benthic_12s","surface_16s", "benthic_16s")
n<-1
for(i in unique(subDat$lake)[-18]){
  print(n)
  n<-n+1
  print(i)
#  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
  par(mfrow=c(4,2))
  z<-0
for(t in geartotest){  
  print(t)
allDat2<-subDat%>%filter(gear == t)
  

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

Data_env<-as.character(Datanew$gear)
Data_env<-Data_env[rowstokeep]
Data_env<-as.factor(Data_env)
Sp_Dat<- data.frame(Datasp,data.xy)
sp_lm<-gls()



geartotest<-c("surface_12s", "benthic_12s","surface_16s", "benthic_16s")
n<-1
for(i in unique(subDat$lake)[-18]){
  print(n)
  n<-n+1
  print(i)
  #  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
#  par(mfrow=c(4,2))
  z<-0
  for(t in geartotest){  
    print(t)
    allDat2<-subDat%>%filter(gear == t)
    
    
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
    data.xy<-data.xy[rowstokeep,]
    
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
    
    
    # break
  }
  break
  mtext(i, side = 3, line = -2, outer = TRUE)
  #  dev.off()
}

xy_dist
transition(lakes)
shortestPath(lakes,)
lakes@data[,1]
r<-raster(ncol=250, nrow=250)
extent(r)<-extent(lakes)
rp<-rasterize(lakes,r,field=1)
?rasterize
plot(rp)

r<-system.file("Lake_Polygons.shp",package = "raster")
r
str(AustinPoly)

library(raster)
library(sf)
# lat
AustinPoly<-sf::st_read("Lake_Austin2.kml")
# AustinPoly
# grid<-st_make_grid(AustinPoly,cellsize = 0.0005, square=T) # here I specify the cell to be of 100x100 meters 
# grid
# plot(grid)
# #AustinPoly<-st_make_valid(AustinPoly)
# grid_int<-st_intersection(grid,AustinPoly)
# grid_sp<-as(grid_int, "Spatial")
# r<-raster(grid_sp)
# res(r)<-c(.0001,.0001)
# r
# plot(grid_int)
# str(grid_int)



#install.packages("fasterize")

library(sp)
library(sf)
files<-list.files(pattern = "*.kml")
for(z in 1:22){
  
  rm(AustinPoly_raster)
  
  
#AustinPoly<-sf::st_read("Lake_Austin2.kml")
  if(file.size(files[z])>20000){
  rv<-2
    } else {
    rv<-1
  }
  print(z)
  t<-unique(subDat$lake)[z]
  print(t)
  #  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
  
  Datanew<-subDat%>%filter(lake == t)
  data.xy<-Datanew[,4:5]
  zone<-floor((data.xy[1,2] + 180) / 6) + 1
  print(zone)
  AustinPoly<-sf::st_read(files[z])
#  str(AustinPoly)
#  st_coordinates(AustinPoly)[,1]
  
#str(AustinPoly)
AustinPoly5<-st_transform(AustinPoly,crs=paste0("+proj=utm +zone=",zone," +datum=WGS84"))
xymat<-as.matrix(data.xy[,c(2,1)])
xy.utm<-rgdal::project(xymat, paste0("+proj=utm +zone=",zone," +datum=WGS84")) 
sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)

xy.utm<-xy.utm[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,]
sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)
AustinPoly5<-st_zm(AustinPoly5)
template <- raster::raster(AustinPoly5, res = rv) #100m resolution
AustinPoly_raster <- fasterize(AustinPoly5, template)
#plot(AustinPoly_raster)





plot(AustinPoly_raster)
points(xy.utm,pch=3)

saveRDS(AustinPoly_raster,paste0(getwd(),substr(files[z],1,9),"r.RDS"))
saveRDS(xy.utm,paste0(getwd(),substr(files[z],1,9),"points.RDS"))

distr<-lcp_over_surface(xy.utm,xy.utm,AustinPoly_raster)
saveRDS(distr,paste0(nam,"dist2.RDS"))
#library(sp)
}
#data.xy
#?projectExtent
#project(AustinPoly,"+proj=utm +zone=16 +datum=WGS84")

#str(AustinPoly5)
AustinPoly5<-st_zm(AustinPoly5)
template <- raster(AustinPoly5, res = 25) #100m resolution
library(fasterize)
AustinPoly_raster <- fasterize(AustinPoly5, template)
plot(AustinPoly_raster)
points(xy.utm,pch=3)


xymat<-as.matrix(data.xy[,c(2,1)])
xy.utm<-project(xymat, "+proj=utm +zone=16 +datum=WGS84") 
pkgload:::unload("dplyr")

library(flapper)
library(terra)
library(raster)
plot(AustinPoly_raster)
library(flapper)
lcp_over_surface(xy.utm,xy.utm,AustinPoly_raster,goal=1)
str(AustinPoly_raster)
AustinPoly_raster@srs<-AustinPoly_raster@crs

template2 <- projectExtent(AustinPoly_raster, crs=("+proj=utm +zone=16 +datum=WGS84")) #100m resolution
resol
AustinPoly_raster2<- projectRaster(AustinPoly_raster,template2,100,crs=("+proj=utm +zone=16 +datum=WGS84") )
projectExtent(object, crs)
plot(AustinPoly_raster2)
templaten <- raster(AustinPoly_raster2, res = 100) #100m resolution

raster::resample(r)
r<-AustinPoly_raster2

lcp_over_surface(xy.utm,xy.utm,AustinPoly_raster,goal=1)

#setwd("E:")

saveRDS(AustinPoly_raster,"r.RDS")
saveRDS(xy.utm,"xy.RDS")
xymat<-as.matrix(data.xy[,c(2,1)])
xy.utm<-project(xymat, "+proj=utm +zone=16 +datum=WGS84") 


#shortestPath(AustinPoly, my.sf.point[1,], my.sf.point[2,], output = "SpatialLines")

#AustinPoly_raster2<-AustinPoly_raster
#AustinPoly_raster[is.na(AustinPoly_raster)] 
r<-AustinPoly_raster
plot(AustinPoly_raster)
r@data@values



my.sf.point <- st_as_sf(x = data.xy, 
                        coords = c("long", "lat"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
plot(my.sf.point, add=T)
points(data.xy[,c(2,1)],pch=3,col="black")

template <- raster(my.sf.point, res = .00005) #100m resolution
#AustinPoly_raster3 <- fasterize(my.sf.point, template)
r <- gdistance::transition(r, mean, directions = 2,symm=F)
class(r)
str(r)

distance <- shortestPath(r, my.sf.point[1,], my.sf.point[2,], output = "SpatialLines")

remotes::install_github("edwardlavender/flapper")
library(flapper)

xymat<-as.matrix(data.xy[,c(2,1)])
range(xymat[,1])

lcp_over_surface(xymat,xymat,r,goal=1)

for(i in unique(subDat$lake)[-18]){
#  print(n)
#  n<-n+1
  print(i)
  #  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
 
  
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
    print(data.xy[1,])
  #  data.xy<-data.xy[rowstokeep,]
}

for(i in unique(subDat$lake)){
  print(i)
  Datanew<-allDat2%>%filter(lake == i)
  data.xy<-Datanew[,4:5]
  print(data.xy[1,])
}



fruits<-round(runif(100,0,238))
fruits

rain<-15+.1566667*fruits+rnorm(100,0,6)

plot(fruits~rain)
