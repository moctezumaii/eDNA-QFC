library(sp)
library(sf)
library(flapper)
library(fasterize)
library(stars)
library(tidyverse)
setwd("D:/rasters/")
files<-list.files(pattern = "*.kml")
for(z in 1:22){
  
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
  
  saveRDS(AustinPoly_raster,paste0(getwd(),"/",substr(files[z],1,9),"r.RDS"))
  saveRDS(xy.utm,paste0(getwd(),"/",substr(files[z],1,9),"points.RDS"))
  xy.utm<-as.matrix(xy.utm)
  distr<-lcp_over_surface(as.matrix(xy.utm),as.matrix(xy.utm),AustinPoly_raster,combination="matrix")
  saveRDS(distr,paste0(getwd(),"/",substr(files[z],1,9),"dist2.RDS"))
  #library(sp)
}