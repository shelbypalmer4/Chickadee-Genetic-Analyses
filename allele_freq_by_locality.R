#### Generate allele frequency by sampling locality ####

setwd("C:/Users/Shelby Palmer/Desktop/CHICKADEES")
PoecileAlleles <- read.csv("PCR-RFLP_data_genetics_24apr2023.csv")

PoecileAlleleFreq <- data.frame(Locality=unique(PoecileAlleles$Pop))

for(i in 1:length(PoecileAlleleFreq$Locality)) {
  PoecileAlleleFreq$BC_allele[i] <- length(which(PoecileAlleles$Pop==PoecileAlleleFreq$Locality[i] & PoecileAlleles[]==1))
  PoecileAlleleFreq$CA_allele[i] <- length(which(PoecileAlleles$Pop==PoecileAlleleFreq$Locality[i] & PoecileAlleles[]==2))
  PoecileAlleleFreq$BC_allele_freq[i] <- PoecileAlleleFreq$BC_allele[i]/sum(PoecileAlleleFreq$BC_allele[i], PoecileAlleleFreq$CA_allele[i])
  PoecileAlleleFreq$CA_allele_freq[i] <- PoecileAlleleFreq$CA_allele[i]/sum(PoecileAlleleFreq$BC_allele[i], PoecileAlleleFreq$CA_allele[i])
}

View(PoecileAlleleFreq)


