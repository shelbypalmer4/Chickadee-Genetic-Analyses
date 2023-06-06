#### Generate allele frequency by sampling locality ####

setwd(genWDlaptop)
PoecileAlleles <- read.csv("PCR-RFLP_data_genetics_24apr2023.csv")

PoecileAlleleFreq <- data.frame(Locality=unique(PoecileAlleles$Pop))

for(i in 1:length(PoecileAlleleFreq$Locality)) {
  PoecileAlleleFreq$BC_allele[i] <- length(which(PoecileAlleles$Pop==PoecileAlleleFreq$Locality[i] & PoecileAlleles[]==1))
  PoecileAlleleFreq$CA_allele[i] <- length(which(PoecileAlleles$Pop==PoecileAlleleFreq$Locality[i] & PoecileAlleles[]==2))
  PoecileAlleleFreq$BC_allele_freq[i] <- PoecileAlleleFreq$BC_allele[i]/sum(PoecileAlleleFreq$BC_allele[i], PoecileAlleleFreq$CA_allele[i])
  PoecileAlleleFreq$CA_allele_freq[i] <- PoecileAlleleFreq$CA_allele[i]/sum(PoecileAlleleFreq$BC_allele[i], PoecileAlleleFreq$CA_allele[i])
}

View(PoecileAlleleFreq)

# 05 June 2023
# create figure for allele frequencies by locality
list.files()
PoecileLoc <- read.csv(list.files()[13])
View(PoecileLoc)
View(PoecileAlleles)

PoecileLoc <- PoecileLoc[order(PoecileLoc$ind_ID),]
PoecileAlleles <- PoecileAlleles[order(PoecileAlleles$Ind),]

PoecileAlleleLoc <- merge(PoecileLoc, 
                          PoecileAlleles, 
                          by.x = "ind_ID",
                          by.y = "Ind")
PoecileAlleleLoc <- PoecileAlleleLoc[order(PoecileAlleleLoc$pop_ID),]

PoecileAlleleFreq2 <- data.frame(Locality=unique(PoecileAlleleLoc$pop_ID))
for(i in 1:length(PoecileAlleleFreq2$Locality)) {
  PoecileAlleleFreq2$BC_allele_count[i] <- length(which(PoecileAlleleLoc[which(PoecileAlleleLoc$pop_ID==i),][]==1))
  PoecileAlleleFreq2$CA_allele_count[i] <- length(which(PoecileAlleleLoc[which(PoecileAlleleLoc$pop_ID==i),][]==2))
  PoecileAlleleFreq2$BC_allele_freq[i] <- PoecileAlleleFreq2$BC_allele_count[i]/sum(PoecileAlleleFreq2$BC_allele_count[i],
                                                                                    PoecileAlleleFreq2$CA_allele_count[i])
  PoecileAlleleFreq2$CA_allele_freq[i] <- PoecileAlleleFreq2$CA_allele_count[i]/sum(PoecileAlleleFreq2$BC_allele_count[i],
                                                                                    PoecileAlleleFreq2$CA_allele_count[i])
  PoecileAlleleFreq2$Ind_count[i] <- length(which(PoecileAlleleLoc$pop_ID==i))
}

# add coordinates for mapping; use median lat/long per locality
for (i in 1:length(PoecileAlleleFreq2$Locality)) {
  PoecileAlleleFreq2$Latitude[i] <- median(PoecileAlleleLoc$latitude[which(PoecileAlleleLoc$pop_ID==i)])
  PoecileAlleleFreq2$Longitude[i] <- median(PoecileAlleleLoc$longitude[which(PoecileAlleleLoc$pop_ID==i)])
}

# add BC, HZ, CA to specify groups
Loc_Ancestry <- c(rep("BCCH",5),
                  rep("Hybrid Zone",2),
                  rep("CACH",3))
PoecileAlleleFreq2 <- cbind(PoecileAlleleFreq2, Loc_Ancestry)

# making the maps

library(ggplot2)
library(ggmap)

# # getting the map
# register_google(key = "AIzaSyBiYXCZirVz62b5oc9hSkx7CGvzc5p8hVM")
# gmap <- get_map(location = c(min(PoecileAlleleFreq2$Longitude)-0.2,
#                              min(PoecileAlleleFreq2$Latitude)-0.15,
#                              max(PoecileAlleleFreq2$Longitude)+0.2,
#                              max(PoecileAlleleFreq2$Latitude)+0.2), 
#                 zoom = 8,
#                 maptype = "satellite")
gmap <- get_map(location = c(min(PoecileAlleleFreq2$Longitude)-0.2,
                             min(PoecileAlleleFreq2$Latitude)-0.15,
                             max(PoecileAlleleFreq2$Longitude)+0.2,
                             max(PoecileAlleleFreq2$Latitude)+0.2), 
                source = "stamen",
                maptype = "terrain")

# plotting the map with locality data and a line approximating the hybrid zone
palette = c("dodgerblue2", "indianred2", "mediumpurple")
png(filename = "Sampling_Locality_Map.png",
    width=6,
    height=3.25,
    units="in",
    res=1200)
ggmap(gmap) +
  geom_point(data = PoecileAlleleFreq2, 
             aes(x = Longitude, 
                 y = Latitude,
                 col = Loc_Ancestry),
             size = 4,
             alpha = 0.8) +
  geom_text(data = PoecileAlleleFreq2,
              aes(x = Longitude,
                  y = Latitude,
                  label = Ind_count),
            size = 2) +
  geom_curve(aes(x = -96.8,
                 y = 37.4,
                 xend = -94.7,
                 yend = 37.7),
             curvature = 0,
             linewidth = 3,
             alpha = 0.1) +
  geom_curve(aes(x = -94.7,
                 y = 37.7,
                 xend = -93.7,
                 yend = 38.4),
             curvature = 0.25,
             angle = 45,
             linewidth = 3,
             alpha = 0.1) +
  geom_curve(aes(x = -93.7,
                 y = 38.4,
                 xend = -93.05,
                 yend = 38.55),
             curvature = -0.35,
             angle = 135,
             linewidth = 3,
             alpha = 0.1) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_manual(values = palette) +
  guides(color = guide_legend(title = "Species Range")) +
  coord_quickmap()
dev.off()


