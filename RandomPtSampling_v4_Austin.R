#Script to randomly sample points within a lake polygon for eDNA sampling - JDR 7/5/18

#Lake polygons downloaded from:
#https://gis-michigan.opendata.arcgis.com/datasets/lake-polygons
#Lake contours downloaded from:
#https://gis-michigan.opendata.arcgis.com/datasets/65085a22e09e432c92d722d85f467f7a_4

#Download the shapefile in zipped format and extract
#Set working directory to the extracted folder
directory<-getwd()
setwd("Lake_Polygons")

#Several libraries are needed for manipulation of shapefiles and plotting below
library(sp)
library(maptools)
library(maps)
library(prevR)
library(plotrix)
library(TSP)
library(RgoogleMaps)
library(rgdal)

#######################################################
#Change this part of the script for other lakes!
#######################################################
lakename = "Austin Lake"
lakeabbrev = "AUS"
nsamples = 40
lakept = c(-85.551679, 42.176194)   #Center point of the lake from LAGOS database
launchpt = c(-85.556282, 42.180116)   #GPS coordinates of the boat launch
setpts = matrix(data = c(
	-85.561704, 42.181957,
	-85.543352, 42.184364,
	-85.540284, 42.166670), ncol = 2, byrow = TRUE)
moresetpts = NULL
#######################################################


#Read in the Lake polygons
setwd("../Lake_Polygons")

#lakes = readShapePoly("Lake_Polygons")
lakes = readOGR(dsn=paste0(getwd(),"/Lake_Polygons"), layer = "Lake_Polygons")
setwd("../Inland_Lake_Contours")

#depths = readShapeLines("Inland_Lake_Contours")
depths = readOGR(dsn=paste0(getwd(),"/Inland_Lake_Contours"), layer = "Inland_Lake_Contours")

#Which lakes in the dataset match the name we're interested in
whichone = which(lakes$LAKE_NAME == lakename)

#Which of these contain the center point from LAGOS (lakept above l.24)
for(x in 1:length(whichone)) {
	tmp = point.in.SpatialPolygons(lakept[1], lakept[2], lakes[whichone[x],])
	if(tmp == TRUE) {
		thisone = whichone[x]
	}
}

thisone

#Replace the full set of lakes with just the one we're interested in (saving space?)
lakes = lakes[thisone,]

#Figure out which contours apply to this lake
is.in.lake = c()
for(x in 1:length(depths@lines)) {
	is.in.lake[x] = point.in.SpatialPolygons(depths@lines[[x]]@Lines[[1]]@coords[1,1], depths@lines[[x]]@Lines[[1]]@coords[1,2], lakes)
}
lakecontours = which(is.in.lake == TRUE)



#Randomly sample from this polygon nsamples times and convert this object to a data frame
pts = spsample(lakes, nsamples-(length(setpts[,1])+length(moresetpts[,1])), type = "random")
pts = as.data.frame(pts)
pts = rbind(launchpt, pts)
for(r in 1:length(setpts[,1])) {
	pts = rbind(pts,setpts[r,])
}
if(!is.null(moresetpts)) {
	for(r in 1:length(moresetpts[,1])) {
		pts = rbind(pts,moresetpts[r,])
	}
}

#Use functions from the TSP package to compute a least cost route through all sample locations
#This is a version of the "Traveling Salesman Problem"
tsp = TSP(dist(pts))
tour = solve_TSP(tsp, method = "nearest_insertion", start = 1L)

#No longer want an asymmetric route -- this is a simple TSP 
#tsp = insert_dummy(tsp, label = "cut")
#tour = solve_TSP(tsp, method = "2-opt", control = list(rep=10))
#path.tsp = unname(cut_tour(tour,"cut"))

#Create a new object for plotting that has sample IDs and GPS coords
#Randomly assign ~30% of the sample as benthic
type = rep("S", nsamples)
nbenthic = floor(0.3*nsamples)  #30% of samples should be benthic
benth = sample(c(1:nsamples), nbenthic, replace = FALSE)
type[benth] = "B"
#tabl = data.frame(name = c("launch",c(1:nsamples)), desc = c(NA,type), latitude = round(pts[path.tsp,2],6), longitude = round(pts[path.tsp,1],6))
tabl = data.frame(name = c("launch",paste0(lakeabbrev,"-",c(1:nsamples))), desc = c(NA,type), latitude = round(pts[tour,2],6), longitude = round(pts[tour,1],6))

#Write this to a csv file
setwd("..")
write.csv(tabl, row.names = F, quote = F, file = paste0(lakename, "2018.csv"))

#Make a pdf plot for the field
pdf(paste0(lakename,"2", ".pdf"), width = 24, height = 12)
par(mfrow = c(1,2))

#OLD PLOTTING DOWN HERE
plot(lakes, main = lakename, lwd = 2, add = TRUE)
lines(tabl$longitude, tabl$latitude, pch = 1, type = "b", lwd = 0.5, cex = 4)
text(tabl$longitude, tabl$latitude, labels = tabl$name)

#New version uses RgoogleMaps
library(RgoogleMaps)
bb = qbbox(tabl$latitude, tabl$longitude)
tmp = paste0(tabl$latitude, ",", tabl$longitude)
mypath = paste0("|", tmp)
mypath = paste(mypath, collapse = "", sep = "")
mypath = paste0("&path=color:grey|weight:1", mypath)
MyMap = GetMap.bbox(bb$lonR, bb$latR, destfile = paste0(lakename,".png"), path = mypath, GRAYSCALE = TRUE)
PlotOnStaticMap(MyMap, lat = tabl$latitude, lon=tabl$longitude, FUN = points, pch= 1, cex = 4)
PlotOnStaticMap(MyMap, lat = tabl$latitude, lon = tabl$longitude, FUN=text, labels = tabl$name, font = 2, cex = 0.5, add = TRUE)
for(x in lakecontours) {
	PlotOnStaticMap(MyMap, lon = depths@lines[[x]]@Lines[[1]]@coords[,1], lat = depths@lines[[x]]@Lines[[1]]@coords[,2], FUN=lines, lwd = 0.25, add = TRUE)
	PlotOnStaticMap(MyMap, lon = depths@lines[[x]]@Lines[[1]]@coords[1,1], lat = depths@lines[[x]]@Lines[[1]]@coords[1,2], FUN= text, font = 2, cex = 0.25, labels = depths@data$DEPTH[x], add = TRUE)
}
plot.new()
addtable2plot(0,0, tabl, display.rownames = F, display.colnames = F, bty = "o", hlines = T, vlines = T, xpad = 2, ypad = 1)
dev.off()






