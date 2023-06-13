#originally created on: October 23th, 2017
#by: Nick Sard

#ABOUT: This script was written to calculate measures of alpha diversity for the eDNA community matrix

#loading necessary libraries
library(tidyverse)
library(permute)
library(lattice)
library(vegan)

#setting working directory and reading in data


#listing files
#list.files("Input/")

#reading in the community matrix for ST
#df <- read.table("Input/fish.mean.max.clean.community.matrix_genus_no_mock.txt", header = T,sep="\t",stringsAsFactors = F)
df <- read.csv("data/12S_100m_22lakes_eDNA_trad_combined_only_spp_level_110719.csv",header = T)
head(df)
colnames(df)
#df <- df %>% select(lake_name,sample,gear,Alosa_pseudoharengus :Carassius_auratus)
names(df)[1] <- "lake"
names(df)[3] <- "gear" #this is named "gear" coz in traditional data there are different gear types
head(df)

#for each lake calculating the Shannon's, Simpson's, Inverse Shannon's, species richness
my.lakes <- unique(df$lake)
i <- 1
i <- NULL
outH <- NULL
outS <- NULL
out_count <- NULL
for(i in 1:length(my.lakes)){
  
  #getting each lake separately
  df1 <- df %>% filter(lake == my.lakes[i]) %>% select(-lake:-gear)
  df2 <- df %>% filter(lake == my.lakes[i]) %>% select(lake:gear)
  head(df1)
  head(df2)
  
  #keeping only the columns with values > 0
  mycols <- colnames(df1)[colSums(df1)!=0]
  df1 <- df1[,mycols]
  head(df1)
  
  #getting counts of species per lake
  all.spp <- ncol(df1)
  counts1 <- data.frame(lake = my.lakes[i], overall = all.spp,stringsAsFactors = F)
  out_count <- rbind(out_count,counts1)
  
  #now for the stats
  #calculating alpha
  H1 <- diversity(df1, index = "shannon")
  H2 <- diversity(df1, index = "simpson")
  H3 <- diversity(df1, index = "invsimpson")
  H <- as.data.frame(t(rbind(H1,H2,H3)))
#  H$effnum <- df2$gear
  H$effnum <- unlist(row.names(H))
  H$type <- "overall"
  H <- cbind(H,df2)
  H <- gather(H,"alpha","value",1:3)
  H$alpha[H$alpha == "H1"] <- "Shannon"
  H$alpha[H$alpha == "H2"] <- "Simpson"
  H$alpha[H$alpha == "H3"] <- "Inverse Simpson"
  head(H)

  #saving for later
  outH <- rbind(outH,H)
  
  #spp richness
  Sover <- as.data.frame(specnumber(df1))
#  Sover$effnum <- df2$gear
  Sover$effnum <- unlist(row.names(Sover))
  colnames(Sover) <- c("R", "effnum")
  Sover <- Sover[,c(2,1)]
  Sover$type <- "overall"
  Sover <- cbind(Sover,df2)
  head(Sover)
  
  #saving for later
  outS <- rbind(outS,Sover)
  
}
head(outH)
head(outS)
head(out_count)
head(df)

#writing these to file
write.table(x = outH, file = "Output/12S_0diff_eDNA_alpha.diversity.per.effort.per.lake.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.csv(outH, "Output/12S_0diff_eDNA_alpha.diversity.per.effort.per.lake.csv", quote = FALSE, row.names = FALSE)
write.table(x = outS, file = "Output/12S_0diff_eDNA_spp.per.effort.per.lake.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.csv(outS, "Output/12S_0diff_eDNA_spp.per.effort.per.lake.csv", quote = FALSE, row.names = FALSE)
write.table(x = out_count, file = "Output/12S_0diff_eDNA_spp.counts.per.lake.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.csv(out_count, "Output/12S_0diff_eDNA_spp.counts.per.lake.csv", quote = FALSE, row.names = FALSE)

#fin!
