#### Calculating LD for the loci from McQuillan et al. 2017 ####

# I get a warning saying the package is now obsolete, but presumably this just means it is no longer being updated but still is functional
library(genetics)
setwd("C:/Users/Shelby Palmer/Desktop/CHICKADEES")
Poecile_Gen<-read.csv("PCR-RFLP_data_for_poppr_31mar2023.csv")
View(Poecile_Gen)
# remove the genalex bullshit
PoecileGen2 <- data.frame(Poecile_Gen[3:length(Poecile_Gen$X9),])
colnames(PoecileGen2) <- c(as.character(Poecile_Gen[2,]))

# getting locus names from the df
loci <- c(colnames(PoecileGen2[,which(1:ncol(PoecileGen2)%%2==1)]))
loci <- loci[2:length(loci)]
# making duplicate names for allele 2
loci2 <- c()
for (i in 1:length(loci)) {
  loci2[i] <- paste(loci[i], "2", sep = "_")
}

colnames(PoecileGen2)[seq(4, ncol(PoecileGen2), 2)] <- loci2

write.csv(PoecileGen2, "PCR-RFLP_data_genetics_24apr2023.csv")
  
# make each locus into a genotype object (should be able to do this iteratively but I gave up on myself)
c0p183_gen <- genotype(c(paste(PoecileGen2$c0p183,
                      PoecileGen2$c0p183_2,
                      sep = "/")))
c0p184_gen <- genotype(c(paste(PoecileGen2$c0p184,
                      PoecileGen2$c0p184_2,
                      sep = "/")))
c0p238_gen <- genotype(c(paste(PoecileGen2$c0p238,
                      PoecileGen2$c0p238_2,
                      sep = "/")))
c0p251_gen <- genotype(c(paste(PoecileGen2$c0p251,
                      PoecileGen2$c0p251_2,
                      sep = "/")))
c0p283_gen <- genotype(c(paste(PoecileGen2$c0p283,
                      PoecileGen2$c0p283_2,
                      sep = "/")))
c0p303_gen <- genotype(c(paste(PoecileGen2$c0p303,
                      PoecileGen2$c0p303_2,
                      sep = "/")))
c0p356_gen <- genotype(c(paste(PoecileGen2$c0p356,
                      PoecileGen2$c0p356_2,
                      sep = "/")))
c0p373_gen <- genotype(c(paste(PoecileGen2$c0p373,
                      PoecileGen2$c0p373_2,
                      sep = "/")))
c0p628_gen <- genotype(c(paste(PoecileGen2$c0p628,
                      PoecileGen2$c0p628_2,
                      sep = "/")))

c0pGens <- data.frame(c0p183_gen, c0p184_gen, c0p238_gen, c0p251_gen, c0p283_gen, c0p303_gen, c0p356_gen, c0p373_gen, c0p628_gen)

c0pGens[c0pGens == "-1/-1"] <- NA
c0pGens <- droplevels(c0pGens)
c0pGens <- makeGenotypes(c0pGens)

# all individuals
c0pLD <- LD(c0pGens)

c0pLD$D
c0pLD$`D'`

# subsetting to only HZ genotypes
c0pHZGens <- c0pGens[1:34,]

c0pHZLD <- LD(c0pHZGens)

c0pHZLD$D
c0pHZLD$`D'`

# Do marker pairs on the same chromosome have higher levels of LD than those on different chromosomes?

## D
# 183, 238, 303
Chr1Pairs <- c(c0pHZLD$D[1,3], 
               c0pHZLD$D[1,6], 
               c0pHZLD$D[3,6])
# 184, 373
Chr3Pairs <- c0pHZLD$D[2,8]
# 283, 356, 628
Chr21Pairs <- c(c0pHZLD$D[5,7],
                c0pHZLD$D[5,9],
                c0pHZLD$D[7,9])

SameChrPairs_D <- c(Chr1Pairs, Chr3Pairs, Chr21Pairs)
DiffChrPairs_D <- c0pHZLD$D[!(c0pHZLD$D %in% SameChrPairs_D)]
DiffChrPairs_D <- DiffChrPairs_D[!is.na(DiffChrPairs_D)]

## D'
# 183, 238, 303
Chr1Pairs <- c(c0pHZLD$`D'`[1,3], 
               c0pHZLD$`D'`[1,6], 
               c0pHZLD$`D'`[3,6])
# 184, 373
Chr3Pairs <- c0pHZLD$`D'`[2,8]
# 283, 356, 628
Chr21Pairs <- c(c0pHZLD$`D'`[5,7],
                c0pHZLD$`D'`[5,9],
                c0pHZLD$`D'`[7,9])

SameChrPairsDPrime <- c(Chr1Pairs, Chr3Pairs, Chr21Pairs)
DiffChrPairsDPrime <- c0pHZLD$`D'`[!(c0pHZLD$`D'` %in% SameChrPairsDPrime)]
DiffChrPairsDPrime <- DiffChrPairsDPrime[!is.na(DiffChrPairsDPrime)]



## r
# 183, 238, 303
Chr1Pairs <- c(c0pHZLD$r[1,3], 
               c0pHZLD$r[1,6], 
               c0pHZLD$r[3,6])
# 184, 373
Chr3Pairs <- c0pHZLD$r[2,8]
# 283, 356, 628
Chr21Pairs <- c(c0pHZLD$r[5,7],
                c0pHZLD$r[5,9],
                c0pHZLD$r[7,9])

SameChrPairs_r <- c(Chr1Pairs, Chr3Pairs, Chr21Pairs)
DiffChrPairs_r <- c0pHZLD$r[!(c0pHZLD$r %in% SameChrPairs_r)]
DiffChrPairs_r <- DiffChrPairs_r[!is.na(DiffChrPairs_r)]


# make a dataframe with values for each available calculation of LD for same- and different-chromosome pairs

HZLDscores <- data.frame(Same_Chr = c(rep("Y", 7),
                                      rep("N", 29)),
                         D = c(SameChrPairs_D, DiffChrPairs_D),
                         D_Prime = c(SameChrPairsDPrime, DiffChrPairsDPrime),
                         r = c(SameChrPairs_r, DiffChrPairs_r))
par(mar=c(5,5,2,5))
stripchart(HZLDscores$D~HZLDscores$Same_Chr, 
           vertical=T,
           pch=16,
           method="jitter",
           jitter=0.01,
           xlab="Same Chromosome?",
           ylab="D",
           col=c("blue", "red"))

stripchart(HZLDscores$D_Prime~HZLDscores$Same_Chr, 
           vertical=T)

stripchart(HZLDscores$r~HZLDscores$Same_Chr, 
           vertical=T)

# base R plotting fucking sucks
library(ggplot2)
library(cowplot)

# D
ggplot(HZLDscores, 
       aes(x=Same_Chr, 
                       y=D, 
                       color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=2) +
  theme(legend.position="none") +
  theme_cowplot(12)

# D'
ggplot(HZLDscores, 
       aes(x=Same_Chr, 
           y=D_Prime, 
           color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=2) +
  theme(legend.position="none") +
  theme_cowplot(12)

# r
ggplot(HZLDscores, 
       aes(x=Same_Chr, 
           y=r, 
           color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=2) +
  theme(legend.position="none") +
  theme_cowplot(12)

