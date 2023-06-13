#Load required packages

library(vegan)
library(ggplot2)
library(viridis)
library(cowplot)
library(gridExtra)

############################################
# Whole community alpha diversity analysis #
############################################

df1<-read.csv("plant.community.lake.csv", header = TRUE)

# Shannon's H'
H1 <- diversity(df1)

# Observed Richness
richness1 <- specnumber(df1)  



############################################
#       Terrestrial alpha diversity        #
############################################

df2<-read.csv("terrestrial.community.lake.csv", header = TRUE)

# Shannon's H'
H2 <- diversity(df2)

# Observed Richness
richness2 <- specnumber(df2)  



############################################
#         Wetland alpha diversity          #
############################################

df3<-read.csv("wetland.community.lake.csv", header = TRUE)

# Shannon's H'
H3 <- diversity(df3)

# Observed Richness
richness3 <- specnumber(df3)  


############################################
#         Wetland alpha diversity          #
############################################

df4<-read.csv("aquatic.community.lake.csv", header = TRUE)

# Shannon's H'
H4 <- diversity(df4)

# Observed Richness
richness4 <- specnumber(df4)  


############################################
#           Algae alpha diversity          #
############################################

df5<-read.csv("algae.community.lake.csv", header = TRUE)

# Shannon's H'
H5 <- diversity(df5)

# Observed Richness
richness5 <- specnumber(df5)  






