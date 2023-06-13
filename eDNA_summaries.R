##summarize eDNA and traditional data
##created by L. Nathan on 12/02/2019

require(gtools)
require(doBy)
require(reshape)
require(ggplot2)

#################################
#####     Read in Data        ###
#################################
#Set working directory and read in data file
#setwd("C:/Users/nathanl/OneDrive - State of Michigan DTMB/eDNA/Metabarcode_occupancy/")

#dat12s<-read.csv("data/12S_100m_22lakes_eDNA_trad_combined_only_spp_level_110719.csv")
#dat16s<-read.csv("data/16S_100m_22lakes_eDNA_trad_combined_only_spp_level_110719.csv")


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
table(subDat$gear)
#add gene to gear label

subDat$gearLab<-NA
subDat$gearLab[grep("_",subDat$sample)]<-paste0(subDat$gear[grep("_",subDat$sample)],"_",subDat$gene[grep("_",subDat$sample)])
subDat$gearLab[which(is.na(subDat$gearLab))]<-as.character(subDat$gear[which(is.na(subDat$gearLab))])
table(subDat$gearLab)
#replace gear label
subDat$gear<-subDat$gearLab
#convert to factor
subDat$gear<-as.factor(subDat$gear)
#add ID for gear
subDat$gearID<-as.integer(subDat$gear)

#reorder
names(subDat)
subDat<-subDat[,c(1:7,116,117,8:115)]

#change to P/A
detDat<-subDat
detDat[,5:117][detDat[,5:117]>0] <- 1

#################################
#####     Summaries        ###
#################################
####detections by lake and species
#melt
meltDat<-melt(detDat[,c(1,10:103)],id.vars = "lake")
#summary of species detections by lake
sumDetbyLake<-summaryBy(value~lake+variable,meltDat,FUN=sum)
#how many lakes was each species detected at?
spByLake<-as.data.frame(table(sumDetbyLake[which(sumDetbyLake$value.sum>0),"variable"]))
table(spByLake$Freq)
#summary of species detections
sumDetbySp<-summaryBy(value~variable,meltDat,FUN=sum)
table(sumDetbySp$value.sum)
hist(sumDetbySp$value.sum)

####detections by sample
detDat$numSp<-rowSums(detDat[,10:103])
detDat$gear<-factor(as.factor(detDat$gear),levels(as.factor(detDat$gear))[c(1,9,2,10,3:8,11)])
ggplot(detDat)+geom_boxplot(aes(x=gear,y=numSp))+theme_classic()

##subset to AIS
AISdat<-detDat[,which(names(detDat)%in%c("lake","gear","Alosa_pseudoharengus",
                                         "Cyprinus_carpio","Carassius_auratus",
                                         "Margariscus_margarita","Morone_americana",
                                         "Mylopharyngodon_piceus","Neogobius_melanostomus",
                                         "Petromyzon_marinus","Osmerus_mordax"
                                         ))]
AISdat$numSp<-rowSums(AISdat[,3:11])
####detections by lake and species
#melt
AISmeltDat<-melt(AISdat[,c(1,3:11)],id.vars = "lake")
#summary of species detections by lake
sumAISbyLake<-summaryBy(value~lake+variable,AISmeltDat,FUN=sum)
castDat<-cast(sumAISbyLake,lake~variable)

#look at alewife data
alewife<-detDat[which(detDat$lake%in%c("Pentwater Lake","Dumont Lake")),c("lake","gear","Alosa_pseudoharengus")]
alewifeCast<-cast(alewife,lake~gear,fun.aggregate = sum)


#################################
##### Proportion of sequences ###
#################################
eDNAdat<-droplevels(allDat[grep("_",allDat$gear),])
table(eDNAdat$gear)
#sum sequences by sample
eDNAdat$totalReads<-rowSums(eDNAdat[,8:length(eDNAdat)],na.rm = T)
#sum AIS sequenes by sample
invasives<-c("Cyprinus","Hypophthalmichthys","Mylopharyngodon","Carassius","Neogobius","Osmerus","Alosa","Ctenopharyngodon")
colMatch<-names(eDNAdat)[grep(paste(invasives,collapse="|"),names(eDNAdat))]
eDNAdat$AISreads<-rowSums(eDNAdat[,colMatch],na.rm = T)
#sum by lake
seqSum<-summaryBy(as.formula(paste(paste(colMatch,collapse="+"),"+totalReads+AISreads~lake")),eDNAdat,FUN=sum)

#remove alewife (only used in traditional gear)
names(seqSum)
seqSum<-seqSum[,-2]
seqSum$propAIS<-round(seqSum$AISreads.sum/seqSum$totalReads.sum,2)
