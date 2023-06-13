library(sp)
library(sf)
library(flapper)
library(fasterize)
library(stars)
library(tidyverse)
library(gtools)
library(vegan)
library(nlme)
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


setwd("D:/rasters/")
files<-list.files(pattern = "*.kml")



for(z in 3:22){
 # par(mfrow=c(4,2))
  rm(AustinPoly_raster)
  
  
  #AustinPoly<-sf::st_read("Lake_Austin2.kml")
  # if(file.size(files[z])>20000){
  #   rv<-3
  # } else {
  #   rv<-1
  # }
  rv<-5
  print(z)
  t<-unique(subDat$lake)[z]
  print(t)
  
  
  #  pdf(paste0(i,"_mantell2.PDF"), 8.5,11)
  
  Datanew<-subDat%>%filter(lake == t)
  Datanew1<-Datanew%>%filter(gear =="surface_12s")
  data.xy<-Datanew1[,4:5]
  Datasr<-Datanew1[,8:ncol(Datanew1)]
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
  
  Datasr<-Datasr[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,] #NEW TRY
  xy.utm<-xy.utm[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,]
  
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)
  AustinPoly5<-st_zm(AustinPoly5)
  template <- raster::raster(AustinPoly5, res = rv) #100m resolution
  AustinPoly_raster <- fasterize(AustinPoly5, template,fun="sum")
  #plot(AustinPoly_raster)
  newRast<-terra::rast(AustinPoly_raster)
  
  
  
  
  plot(AustinPoly_raster)
  points(xy.utm,pch=3,col="black")
  extracted<-terra::extract(AustinPoly_raster,xy.utm)
  if(NA %in% extracted){
    print(paste0(t, " has points outside the raster"))
    orp<-which(is.na(extracted))
    points(xy.utm[orp,],pch=3,col="red",cex=3)
    mt<-matrix(data=c(rep(-5:5,11),rep(-5:5,each=11)),ncol=2)
    for(i in orp){
      #  print(i)
      original<-xy.utm[i,]
      for (new in 1:(nrow(mt))){
        
        if (is.na(terra::extract(AustinPoly_raster,xy.utm)[i])){
          #  print(new)
          #   print(i)
          xy.utm[i,]<-original+mt[new,]
          #    print(xy.utm[i,])
        } else {
          
          break}
      }
      
    }
    extracted<-terra::extract(AustinPoly_raster,xy.utm)
    if(NA %in% extracted){
      print(paste0(t, " STILL has points outside the raster")) }
    # print(extracted)
  } 
  #p2<-rasterToPolygons(AustinPoly_raster)
  #p2<-st_as_sf(p2)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==1)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==0)
  #plot(p2)
  
  print("NUMBER #1")
  
  saveRDS(AustinPoly_raster,paste0(getwd(),"/",substr(files[z],1,9),"r1.RDS"))
  saveRDS(xy.utm,paste0(getwd(),"/",substr(files[z],1,9),"points1.RDS"))
  xy.utm<-as.matrix(xy.utm)
  distr<-lcp_over_surface(as.matrix(xy.utm),as.matrix(xy.utm),AustinPoly_raster,combination="matrix")
  saveRDS(distr,paste0(getwd(),"/",substr(files[z],1,9),"dist_xy_1.RDS"))
  sr_dist<-dist(Datasr)
  saveRDS(sr_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_sr_1.RDS"))
  #library(sp)
  ###
  comm_dist <- vegdist(Datasr,method="manhattan")
  saveRDS(comm_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_comm_1.RDS"))
  
  distlcp<-as.dist(distr$dist_lcp)
  
  comm_mantel1 <- mantel(distlcp, comm_dist)
  sr_corlog1 <- mantel.correlog(sr_dist, distlcp)
  comm_corlog1 <- mantel.correlog(comm_dist, distlcp)
  
  mantelR<-list(comm_mantel1,sr_corlog1,comm_corlog1)
  saveRDS(mantelR,paste0(getwd(),"/",substr(files[z],1,9),"mantel_1.RDS"))
  
  
  #### N2
  print("NUMBER #2")
  Datanew2<-Datanew%>%filter(gear =="benthic_12s")
  data.xy<-Datanew2[,4:5]
  Datasr<-Datanew2[,8:ncol(Datanew1)]
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
  
  Datasr<-Datasr[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,] #NEW TRY
  
  xy.utm<-xy.utm[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,]
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)
  AustinPoly5<-st_zm(AustinPoly5)
  template <- raster::raster(AustinPoly5, res = rv) #100m resolution
  AustinPoly_raster <- fasterize(AustinPoly5, template,fun="sum")
  #plot(AustinPoly_raster)
  newRast<-terra::rast(AustinPoly_raster)
  
  
  
  
  plot(AustinPoly_raster)
  points(xy.utm,pch=3,col="black")
  extracted<-terra::extract(AustinPoly_raster,xy.utm)
  if(NA %in% extracted){
    print(paste0(t, " has points outside the raster"))
    orp<-which(is.na(extracted))
    points(xy.utm[orp,],pch=3,col="red",cex=3)
    mt<-matrix(data=c(rep(-5:5,11),rep(-5:5,each=11)),ncol=2)
    for(i in orp){
      #  print(i)
      original<-xy.utm[i,]
      for (new in 1:(nrow(mt))){
        
        if (is.na(terra::extract(AustinPoly_raster,xy.utm)[i])){
          #  print(new)
          #   print(i)
          xy.utm[i,]<-original+mt[new,]
          #    print(xy.utm[i,])
        } else {
          
          break}
      }
      
    }
    extracted<-terra::extract(AustinPoly_raster,xy.utm)
    if(NA %in% extracted){
      print(paste0(t, " STILL has points outside the raster")) }
    # print(extracted)
  } 
  #p2<-rasterToPolygons(AustinPoly_raster)
  #p2<-st_as_sf(p2)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==1)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==0)
  #plot(p2)
  
  saveRDS(AustinPoly_raster,paste0(getwd(),"/",substr(files[z],1,9),"r2.RDS"))
  saveRDS(xy.utm,paste0(getwd(),"/",substr(files[z],1,9),"points2.RDS"))
  xy.utm<-as.matrix(xy.utm)
  distr<-lcp_over_surface(as.matrix(xy.utm),as.matrix(xy.utm),AustinPoly_raster,combination="matrix")
  saveRDS(distr,paste0(getwd(),"/",substr(files[z],1,9),"dist_xy_2.RDS"))
  sr_dist<-dist(Datasr)
  saveRDS(sr_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_sr_2.RDS"))
  #
  ###
  comm_dist <- vegdist(Datasr,method="manhattan")
  saveRDS(comm_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_comm_2.RDS"))
  
  distlcp<-as.dist(distr$dist_lcp)
  
  comm_mantel2 <- mantel(distlcp, comm_dist)
  sr_corlog2 <- mantel.correlog(sr_dist, distlcp)
  comm_corlog2 <- mantel.correlog(comm_dist, distlcp)
  
  mantelR<-list(comm_mantel2,sr_corlog2,comm_corlog2)
  saveRDS(mantelR,paste0(getwd(),"/",substr(files[z],1,9),"mantel_2.RDS"))
  
  
  
  
  
  #### 3
  Datanew1<-Datanew%>%filter(gear =="surface_16s")
  data.xy<-Datanew1[,4:5]
  Datasr<-Datanew1[,8:ncol(Datanew1)]
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
  
  Datasr<-Datasr[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,] #NEW TRY
  
  xy.utm<-xy.utm[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,]
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)
  AustinPoly5<-st_zm(AustinPoly5)
  template <- raster::raster(AustinPoly5, res = rv) #100m resolution
  AustinPoly_raster <- fasterize(AustinPoly5, template,fun="sum")
  #plot(AustinPoly_raster)
  newRast<-terra::rast(AustinPoly_raster)
  
  
  
  
  plot(AustinPoly_raster)
  points(xy.utm,pch=3,col="black")
  extracted<-terra::extract(AustinPoly_raster,xy.utm)
  if(NA %in% extracted){
    print(paste0(t, " has points outside the raster"))
    orp<-which(is.na(extracted))
    points(xy.utm[orp,],pch=3,col="red",cex=3)
    mt<-matrix(data=c(rep(-5:5,11),rep(-5:5,each=11)),ncol=2)
    for(i in orp){
      #  print(i)
      original<-xy.utm[i,]
      for (new in 1:(nrow(mt))){
        
        if (is.na(terra::extract(AustinPoly_raster,xy.utm)[i])){
          #  print(new)
          #   print(i)
          xy.utm[i,]<-original+mt[new,]
          #    print(xy.utm[i,])
        } else {
          
          break}
      }
      
    }
    extracted<-terra::extract(AustinPoly_raster,xy.utm)
    if(NA %in% extracted){
      print(paste0(t, " STILL has points outside the raster")) }
    # print(extracted)
  } 
  #p2<-rasterToPolygons(AustinPoly_raster)
  #p2<-st_as_sf(p2)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==1)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==0)
  #plot(p2)
  
  saveRDS(AustinPoly_raster,paste0(getwd(),"/",substr(files[z],1,9),"r3.RDS"))
  saveRDS(xy.utm,paste0(getwd(),"/",substr(files[z],1,9),"points3.RDS"))
  xy.utm<-as.matrix(xy.utm)
  distr<-lcp_over_surface(as.matrix(xy.utm),as.matrix(xy.utm),AustinPoly_raster,combination="matrix")
  saveRDS(distr,paste0(getwd(),"/",substr(files[z],1,9),"dist_xy_3.RDS"))
  sr_dist<-dist(Datasr)
  saveRDS(sr_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_sr_3.RDS"))
  #library(sp)
  
  ###
  comm_dist <- vegdist(Datasr,method="manhattan")
  saveRDS(comm_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_comm_3.RDS"))
  
  distlcp<-as.dist(distr$dist_lcp)
  
  comm_mantel3 <- mantel(distlcp, comm_dist)
  sr_corlog3 <- mantel.correlog(sr_dist, distlcp)
  comm_corlog3 <- mantel.correlog(comm_dist, distlcp)
  
  mantelR<-list(comm_mantel3,sr_corlog3,comm_corlog3)
  saveRDS(mantelR,paste0(getwd(),"/",substr(files[z],1,9),"mantel_3.RDS"))
  
  
  
  #### N4
  Datanew4<-Datanew%>%filter(gear =="benthic_16s")
  data.xy<-Datanew4[,4:5]
  Datasr<-Datanew4[,8:ncol(Datanew1)]
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
  
  Datasr<-Datasr[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,] #NEW TRY
  
  xy.utm<-xy.utm[point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1,]
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==1)
  sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(AustinPoly5)[,1], st_coordinates(AustinPoly5)[,2])==0)
  AustinPoly5<-st_zm(AustinPoly5)
  template <- raster::raster(AustinPoly5, res = rv) #100m resolution
  AustinPoly_raster <- fasterize(AustinPoly5, template,fun="sum")
  #plot(AustinPoly_raster)
  newRast<-terra::rast(AustinPoly_raster)
  
  
  
  
  plot(AustinPoly_raster)
  points(xy.utm,pch=3,col="black")
  extracted<-terra::extract(AustinPoly_raster,xy.utm)
  if(NA %in% extracted){
    print(paste0(t, " has points outside the raster"))
    orp<-which(is.na(extracted))
    points(xy.utm[orp,],pch=3,col="red",cex=3)
    mt<-matrix(data=c(rep(-5:5,11),rep(-5:5,each=11)),ncol=2)
    for(i in orp){
      #  print(i)
      original<-xy.utm[i,]
      for (new in 1:(nrow(mt))){
        
        if (is.na(terra::extract(AustinPoly_raster,xy.utm)[i])){
          #  print(new)
          #   print(i)
          xy.utm[i,]<-original+mt[new,]
          #    print(xy.utm[i,])
        } else {
          
          break}
      }
      
    }
    extracted<-terra::extract(AustinPoly_raster,xy.utm)
    if(NA %in% extracted){
      print(paste0(t, " STILL has points outside the raster")) }
    # print(extracted)
  } 
  #p2<-rasterToPolygons(AustinPoly_raster)
  #p2<-st_as_sf(p2)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==1)
  #sum(point.in.polygon(xy.utm[,1],xy.utm[,2], st_coordinates(p2)[,1], st_coordinates(p2)[,2])==0)
  #plot(p2)
  
  saveRDS(AustinPoly_raster,paste0(getwd(),"/",substr(files[z],1,9),"r4.RDS"))
  saveRDS(xy.utm,paste0(getwd(),"/",substr(files[z],1,9),"points4.RDS"))
  xy.utm<-as.matrix(xy.utm)
  distr<-lcp_over_surface(as.matrix(xy.utm),as.matrix(xy.utm),AustinPoly_raster,combination="matrix")
  saveRDS(distr,paste0(getwd(),"/",substr(files[z],1,9),"dist_xy_4.RDS"))
  sr_dist<-dist(Datasr)
  saveRDS(sr_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_sr_4.RDS"))
  
  ###
  comm_dist <- vegdist(Datasr,method="manhattan")
  saveRDS(comm_dist,paste0(getwd(),"/",substr(files[z],1,9),"dist_comm_4.RDS"))
  
  distlcp<-as.dist(distr$dist_lcp)
  
  comm_mantel4 <- mantel(distlcp, comm_dist)
  sr_corlog4 <- mantel.correlog(sr_dist, distlcp)
  comm_corlog4 <- mantel.correlog(comm_dist, distlcp)
  
  mantelR<-list(comm_mantel4,sr_corlog4,comm_corlog4)
  saveRDS(mantelR,paste0(getwd(),"/",substr(files[z],1,9),"mantel_4.RDS"))
  
  pdf(paste0(getwd(),"/",substr(files[z],1,9),"mantel.pdf"),width=5,height=8)
  
  par(mfrow=c(4,2))
  plot(sr_corlog1)
  mtext(side=3, paste0('Species Richness S 12S'))
  plot(comm_corlog1)
  mtext(side=3, paste0('Community Composition S 12S'))
  
  plot(sr_corlog2)
  mtext(side=3, paste0('Species Richness B 12S'))
  plot(comm_corlog2)
  mtext(side=3, paste0('Community Composition B 12S'))
  
  plot(sr_corlog3)
  mtext(side=3, paste0('Species Richness S 16S'))
  plot(comm_corlog3)
  mtext(side=3, paste0('Community Composition S 16S'))
  
  plot(sr_corlog4)
  mtext(side=3, paste0('Species Richness B 16S'))
  plot(comm_corlog4)
  mtext(side=3, paste0('Community Composition B 16S'))
  
  
  mtext(t, side = 3, line = -2, outer = TRUE)
  
  dev.off()
  
}   
