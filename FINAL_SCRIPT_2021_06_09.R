##occupancy model for eDNA metabarcoding data
##Goal: compare detection probability of multiple gear types (eDNA and traditional) across three species
##created by L. Nathan on 4/27/2019, modifed for multiple genes on 9/30/2019

###Assumptions
##based on "Lake-level" closure (i.e. all samples)
##ignores within-lake variation

library(jagsUI)
library(doBy)
library(reshape)
library(ggplot2)
library(tidyr)
library(gtools)

#################################
#####     Read in Data        ###
#################################
#Set working directory and read in data file
#setwd("C:/Users/nathanl/OneDrive - State of Michigan DTMB/eDNA/Pukk_manuscript/Metabarcode_occupancy/")

#save detDat table
#detDat<-read.csv("./data/2019_11_07/detDat_threeSp_byLake_2019-12-03.csv")

#################################
#####    Run Occ Model        ###
#################################

#define model
cat(file = "eDNA_occ_allSp_byLake_2genes.txt","
    model{
    #priors
    for(sp in 1:s){
        psi[sp]~dunif(0,1)
    }
    for(gear in 1:g){
      for(sp in 1:s){
        p[gear,sp]~dunif(0,1)
      }
    }
    

    #likelihood (biological model for true occupancy)
    for(i in 1:n){ #n = number of lakes
      for(sp in 1:s){
        z[i,sp]~dbern(psi[sp])
      }
    }
    #observation model
    for(i in 1:nr){ #nr = total number of samples
      y[i]~dbin(p[gear[i],species[i]],z[lake[i],species[i]])
    }

}")

detDat<-detDat[,-ncol(detDat)]
detDat<-detDat%>%pivot_longer(Alosa_pseudoharengus:Noturus_miurus,names_to="species",values_to = "detect")
zst1<-cast(summaryBy(formula=detect~species+lake,detDat,FUN=max),formula = lake~species,value = "detect.max")
colnames(zst1[-1])

#bundle data
jagsData<-list(n=length(unique(detDat$lake)), #number lakes
               nr=nrow(detDat), #number samples
               g=length(unique(detDat$gear)), #number gears
               s=length(unique(detDat$species)), #number species
               y=detDat$detect, #detection               
               lake=as.integer(as.factor(detDat$lake)), #lake
               gear=as.integer(as.factor(detDat$gear)), #gear
               species=as.integer(as.factor(detDat$species))) #species




#inits function (naive occupancy)
zst<-as.matrix(cast(summaryBy(formula=detect~species+lake,detDat,FUN=max),formula = lake~species,value = "detect.max")[,-1])



inits<-function(){list(z=zst)}

#parameters to estimate
params<-c("psi","p")

#start Gibbs sampler
jagsMod <- jags(data = jagsData,
                inits = inits,
                parameters.to.save = params,
                model.file = "eDNA_occ_allSp_byLake_2genes.txt",
                n.chains = 3,
                n.iter = 5000,
                n.burnin = 1000,
                n.thin = 5)

#################################
#####    Plot Output          ###
#################################

#store posteriors as dataframe
postDF<-as.data.frame(jagsMod$sims.list)


##save detection probabilities
detp<-postDF[,109:(ncol(postDF)-1)]
#create labels
gears<-as.character(levels(detDat$gear))
species<-as.character(levels(as.factor(detDat$species)))
labs<-matrix(nrow=length(gears),ncol=length(species))
for(g in 1:length(gears)){
  for(s in 1:length(species)){
    labs[g,s]<-paste0(gears[g],"-",species[s])
  }
}

names(detp)<-unlist(as.list(labs))

detP<-melt(detp)


detP$variable<-as.character(detP$variable)
names(detp)<-unlist(as.list(labs))
detP<- detP %>% separate(variable,c("gear","species"),sep = "-")

#re order gears
detP$gear<-factor(as.factor(detP$gear),levels(as.factor(detP$gear))[c(1,9,2,10,3:8,11)])


#median and 95% credible intervals
plotDat <- summaryBy(value ~ species+gear,detP,
                     FUN = c(median, function(x) quantile(x, c(.05, .95)),var),
                     fun.names = c("Median", "LCI", "UCI"))
names(plotDat)[3:6]<-c("Median","LCI","UCI","Var")


plotOut<-ggplot(plotDat)+
  geom_point(aes(x=species,y=Median,group=gear,color=gear,shape=gear),
             position = position_dodge(width=.5),size=5)+
  geom_errorbar(aes(x=species,ymin=LCI,ymax=UCI,group=gear,color=gear),
                position = position_dodge(width=.5),width=.25,size=1)+
  ylab("Detection Probability")+xlab("Species")+labs(Color = "Gear")+
  scale_colour_manual(name = "Gear",
                      labels = c("Benthic 12s","Surface 12s","Benthic 16s","Surface 16s",
                                 "Boomshocker","Experimental Gillnet", "Great Lakes Gillnet", "Large Mesh Fyke Net","Seine","Small Mesh Fyke Net",
                                 "Trap Net"),
                      values = c("#FF62BC","#E76BF3","#9590FF","#00B0F6","#00BFC4","#00BF7D","#39B600","#A3A500",
                                 "#D89000","#F8766D","#E10707"))+
  scale_shape_manual(name = "Gear",
                     labels = c("Benthic 12s","Surface 12s","Benthic 16s","Surface 16s",
                                "Boomshocker","Experimental Gillnet", "Great Lakes Gillnet", "Large Mesh Fyke Net","Seine","Small Mesh Fyke Net",
                                "Trap Net"),
                     values = c(19,19,19,19,17,17,17,17,17,17,17))+
  theme_classic()+
  ylim(0,1)+
  coord_flip()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=15),
        legend.position = c(.8,.5),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

plotOut

ggsave(filename="outputs/Figs/occupancyMod_fig_2021_06_08.pdf",
       plot=plotOut,
       device=NULL,
       units="mm",
       dpi=600)

