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
# create figure of parental and HZ sampling locations
list.files()
PoecileLoc <- read.csv("STRUCTURE_data_with_localities.csv")
View(PoecileLoc)
View(PoecileAlleles)

PoecileLoc <- PoecileLoc[order(PoecileLoc$ind_ID),]
PoecileAlleles <- PoecileAlleles[order(PoecileAlleles$Ind),]

PoecileAlleleLoc <- merge(PoecileLoc, 
                          PoecileAlleles, 
                          by.x = "ind_ID",
                          by.y = "Ind")
PoecileAlleleLoc <- PoecileAlleleLoc[order(PoecileAlleleLoc$pop_ID),]
colnames(PoecileAlleleLoc)
PoecileAlleleLoc <- PoecileAlleleLoc[,-c(2:3,11)]
colnames(PoecileAlleleLoc)

PoecileAlleleFreq2 <- data.frame(Locality=unique(PoecileAlleleLoc$pop_ID))
for(i in 1:length(PoecileAlleleFreq2$Locality)) {
  PoecileAlleleFreq2$BC_allele_count[i] <- length(which(PoecileAlleleLoc[which(PoecileAlleleLoc$pop_ID==i),10:length(colnames(PoecileAlleleLoc))][]==1))
  PoecileAlleleFreq2$CA_allele_count[i] <- length(which(PoecileAlleleLoc[which(PoecileAlleleLoc$pop_ID==i),10:length(colnames(PoecileAlleleLoc))][]==2))
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
palette = c("dodgerblue3", "indianred1", "mediumpurple")
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

# 6 June 2023
# Make allele frequency map
# new dataframe with only localities with 5+ samples
PoecileAlleleFreq5plus <- PoecileAlleleFreq2[which(PoecileAlleleFreq2$Ind_count>4),]
# add sampling locality names
Localities <- c("Butler Lake, Bates Co.",
                "Clinton, Henry Co.",
                "Sparrowfoot Park, Henry Co.",
                "Bird Song CA, St. Clair Co.")
PoecileAlleleFreq5plus <- cbind(PoecileAlleleFreq5plus, Localities)

install.packages("scatterpie")
library(ggplot2)
library(ggmap)
library(scatterpie)

# define new map boundaries
hmap <- get_map(location = c(min(PoecileAlleleFreq5plus$Longitude)-0.2,
                             min(PoecileAlleleFreq5plus$Latitude)-0.15,
                             max(PoecileAlleleFreq5plus$Longitude)+0.2,
                             max(PoecileAlleleFreq5plus$Latitude)+0.2), 
                source = "stamen",
                maptype = "terrain")

# make new map
palette = c("dodgerblue3", "indianred1", "mediumpurple")
png(filename = "AlleleFreq_by_Locality_Map.png",
    width=6,
    height=3.25,
    units="in",
    res=1200)
ggmap(hmap) +
  geom_curve(aes(x = -94.51,
                 y = 37.73,
                 xend = -93.53,
                 yend = 38.43),
             curvature = -0.15,
             angle = 45,
             linewidth = 13,
             lineend = "round",
             alpha = 0.1) +
  geom_scatterpie(data = PoecileAlleleFreq5plus, 
                  aes(x = Longitude, 
                      y = Latitude),
                  cols = c("BC_allele_freq",
                           "CA_allele_freq"),
                  pie_scale = 3,
                  color = NA) +
  geom_text(data = PoecileAlleleFreq5plus,
            aes(x = Longitude,
                y = Latitude,
                label = Localities),
            size = 2,
            nudge_x = -0.07,
            nudge_y = -0.05) +
  geom_label(aes(x = -94.5,
                 y = 38.51,
                 label = "Missouri")) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_manual(values = palette,
                    name = "Species Allele",
                    labels = c("BCCH", "CACH")) +
  coord_quickmap() +
  coord_equal()
dev.off()

# map with allele frequencies for all localities. Not going to print yet; just in case Jay asks about this
PoecileAlleleFreq2$pie_scale <- rep(NA)
for (i in 1:length(PoecileAlleleFreq2$pie_scale)) {
  PoecileAlleleFreq2$pie_scale[i] <- PoecileAlleleFreq2$Ind_count[i]*(3/max(PoecileAlleleFreq2$Ind_count))
}
alt_pie_scale <- rep(1:20, 2)
palette = c("dodgerblue3", "indianred1", "mediumpurple")
png(filename = "AlleleFreq_by_Locality_Map_All.png",
    width=6,
    height=3.25,
    units="in",
    res=1200)
ggmap(gmap) +
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
  geom_scatterpie(data = PoecileAlleleFreq2, 
                  aes(x = Longitude, 
                      y = Latitude),
                  cols = c("BC_allele_freq",
                           "CA_allele_freq"),
                  pie_scale = 1.25,
                  color = NA) +
  geom_text(data = PoecileAlleleFreq2,
            aes(x = Longitude,
                y = Latitude,
                label = Ind_count),
            size = 4,
            nudge_x = -0.1) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_manual(values = palette,
                    name = "Species Allele",
                    labels = c("BCCH", "CACH")) +
  coord_quickmap() +
  coord_equal()
dev.off()

# tables for thesis
# add county name column
PoecileAlleleFreq2$cityetc_county <- c("Manhattan, Riley, KS",
                                    "Lawrence, Douglas, KS",
                                    "Lake Quivira, Johnson, KS",
                                    "Leeton, Johnson, MO",
                                    "Butler Lake, Bates, MO",
                                    "Clinton, Henry, MO",
                                    "Sparrowfoot Park, Henry, MO",
                                    "Schell-Osage CA, St. Clair, MO",
                                    "Springfield, Greene, MO",
                                    "Table Rock Lake, Stone, MO")
write.table(cbind(signif(PoecileAlleleFreq2[,4:5],
                         digits = 3),
                  PoecileAlleleFreq2[,6:10]),
            file = "allelefreqbyloc.txt", 
            sep = ",", 
            quote = F,
            row.names = F)

# make new table with KU chickadee loan info
KUPoecile <- PoecileAlleleFreq2[c(1:3,5,8,10),]
write.table(KUPoecile,
            file = "partialKUinfo.txt", 
            sep = ",", 
            quote = F,
            row.names = F)

# summary stats for BCCH and CACH allele freq
BCCHallelefreq <- sum(PoecileAlleleFreq2$CA_allele_count[which(PoecileAlleleFreq2$Locality %in% 1:5)])/sum(PoecileAlleleFreq2$CA_allele_count[which(PoecileAlleleFreq2$Locality %in% 1:5)],
     PoecileAlleleFreq2$BC_allele_count[which(PoecileAlleleFreq2$Locality %in% 1:5)])
BCCHallelefreq

CACHallelefreq <- sum(PoecileAlleleFreq2$CA_allele_count[which(PoecileAlleleFreq2$Locality %in% 8:10)])/sum(PoecileAlleleFreq2$CA_allele_count[which(PoecileAlleleFreq2$Locality %in% 8:10)],
                                                                                                           PoecileAlleleFreq2$BC_allele_count[which(PoecileAlleleFreq2$Locality %in% 8:10)])
BCCHallelefreq

#
#
#### preliminary HZ search map for defense slides ####
searchspots <- data.frame(Location = c("Poague CA", "Clinton City Park", "Sparrowfoot Park", "Cooper Creek ATV Park", "Brush Creek SWMA", "Brownington SWA", "Katy Trail-Calhoun", "Katy Trail-Clinton", "Cole Camp", "Montrose CA", "Truman Lake State Park"),
                          Latitude = c(38.4126, 38.3613, 38.29672, 38.2631, 38.3996, 38.2453, 38.4685, 38.3852, 38.4492, 38.3093, 38.2721),
                          Longitude = c(-93.8384, -93.7922, -93.7342, -93.7466, -93.6480, -93.7346, -93.6233, -93.7572, -93.2011, -93.9712, -93.4431))

imap <- get_map(location = c(min(searchspots$Longitude)-0.05,
                             min(searchspots$Latitude)-0.05,
                             max(searchspots$Longitude)+0.05,
                             max(searchspots$Latitude)+0.05), 
                source = "stamen",
                maptype = "terrain")

# plotting the map with locality data and a line approximating the hybrid zone
setwd(genWDlaptop)
png(filename = "HZsearch_map.png",
    width=6,
    height=3.25,
    units="in",
    res=1200)
ggmap(imap) +
  geom_point(data = searchspots, 
             aes(x = Longitude, 
                 y = Latitude),
             size = 2,
             color = "indianred3") +
  geom_text(data = searchspots,
            aes(x = Longitude,
                y = Latitude,
                label = Location),
            size = 2,
            nudge_x = 0.01,
            nudge_y = -0.008) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_quickmap()
dev.off()