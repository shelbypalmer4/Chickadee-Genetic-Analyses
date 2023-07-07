#### Calculating LD for the loci from McQuillan et al. 2017 ####

# I get a warning saying the package is now obsolete, but presumably this just means it is no longer being updated but still is functional
setwd(genWDlaptop)
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

## r^2
# 183, 238, 303
Chr1Pairs <- c(c0pHZLD$`R^2`[1,3], 
               c0pHZLD$`R^2`[1,6], 
               c0pHZLD$`R^2`[3,6])
# 184, 373
Chr3Pairs <- c0pHZLD$`R^2`[2,8]
# 283, 356, 628
Chr21Pairs <- c(c0pHZLD$`R^2`[5,7],
                c0pHZLD$`R^2`[5,9],
                c0pHZLD$`R^2`[7,9])

SameChrPairs_rsq <- c(Chr1Pairs, Chr3Pairs, Chr21Pairs)
DiffChrPairs_rsq <- c0pHZLD$r[!(c0pHZLD$r %in% SameChrPairs_r)]
DiffChrPairs_rsq <- DiffChrPairs_r[!is.na(DiffChrPairs_r)]

# make a dataframe with values for each available calculation of LD for same- and different-chromosome pairs

HZLDscores <- data.frame(Same_Chr = c(rep("Y", 7),
                                      rep("N", 29)),
                         D = c(SameChrPairs_D, DiffChrPairs_D),
                         D_Prime = c(SameChrPairsDPrime, DiffChrPairsDPrime),
                         r = c(SameChrPairs_r, DiffChrPairs_r),
                         r_sq = c(SameChrPairs_rsq, DiffChrPairs_rsq))
#### 11 May 2023. Let's get means of each group ####

HZLDscores_PL <- HZLDscores[which(HZLDscores$Same_Chr=="Y"),]
HZLDscores_PUL <- HZLDscores[which(HZLDscores$Same_Chr=="N"),]

PL_D_mean <- mean(HZLDscores_PL$D)
PL_Dpr_mean <- mean(HZLDscores_PL$D_Prime)
PL_r_mean <- mean(HZLDscores_PL$r)
PL_rsq_mean <- mean(HZLDscores_PL$r_sq)

PUL_D_mean <- mean(HZLDscores_PUL$D)
PUL_Dpr_mean <- mean(HZLDscores_PUL$D_Prime)
PUL_r_mean <- mean(HZLDscores_PUL$r)
PUL_rsq_mean <- mean(HZLDscores_PUL$r_sq)


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

write.csv(HZLDscores, "HZ_LD_scores.csv")

#### 13 may 2023: LD tests ####
var(HZLDscores$r_sq[which(HZLDscores$Same_Chr=="Y")]) 
var(HZLDscores$r_sq[which(HZLDscores$Same_Chr=="N")]) 
# variances not equal for D or r^2 though...order of magnitude difference. log-transforming helps a little. Welch's it is
ldmodel_D <- t.test(data = HZLDscores,
               D~Same_Chr,
               var.equal = F)
ldmodel_Dpr <- t.test(data = HZLDscores,
                 D_Prime~Same_Chr,
                 var.equal = F)
ldmodel_rsq <- t.test(data = HZLDscores,
                 r_sq~Same_Chr,
                 var.equal = F)
ldmodel_D # p = 0.1492
ldmodel_Dpr # p = 0.1991
ldmodel_rsq # p = 0.6052



# new plot with Welch's t-test values
Dvals <- paste("t(", round(ldmodel_D$parameter,
                           digits = 2), ") = ",
               round(ldmodel_D$statistic,
                          digits = 2),
              "  p = ", round(ldmodel_D$p.value,
                      digits = 3),
               sep = "")
Dprvals <- paste("t(", round(ldmodel_Dpr$parameter,
                             digits = 2), ") = ",
                 round(ldmodel_Dpr$statistic,
                       digits = 2),
                 "  p = ", round(ldmodel_Dpr$p.value,
                                 digits = 3),
                 sep = "")
rsqvals <- paste("t(", round(ldmodel_rsq$parameter,
                             digits = 2), ") = ",
                 round(ldmodel_rsq$statistic,
                       digits = 2),
                 "  p = ", round(ldmodel_rsq$p.value,
                                 digits = 3),
                 sep = "")

library(ggplot2)
library(ggpubr)
theme_set(theme_bw())
# D
Dplot <- ggplot(HZLDscores, 
                aes(x=Same_Chr, 
                    y=D, 
                    color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=3,
              alpha=0.75) +
  scale_x_discrete(labels = c("Physically Unlinked", 
                              "Physically Linked")) +
  geom_segment(aes(x = 0.75,
                   y = PUL_D_mean,
                   xend = 1.25,
                   yend = PUL_D_mean),
               color = "black",
               linewidth = 0.75) +
  geom_segment(aes(x = 1.75,
                   y = PL_D_mean,
                   xend = 2.25,
                   yend = PL_D_mean),
               color = "black",
               linewidth = 0.75) +
  geom_text(aes(x = 2,
                y = max(c0pHZLD$D, na.rm = T),
                label = "c0p356-c0p628"),
            size = 3,
            color = "black",
            nudge_y = -0.006) +
  ggtitle(Dvals) +
  theme(legend.position="none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"))

# D'
DPplot <- ggplot(HZLDscores, 
                 aes(x=Same_Chr, 
                     y=D_Prime, 
                     color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=3,
              alpha=0.75) +
  scale_x_discrete(labels = c("Physically Unlinked", 
                              "Physically Linked")) +
  geom_segment(aes(x = 0.75,
                   y = PUL_Dpr_mean,
                   xend = 1.25,
                   yend = PUL_Dpr_mean),
               color = "black",
               linewidth = 0.75) +
  geom_segment(aes(x = 1.75,
                   y = PL_Dpr_mean,
                   xend = 2.25,
                   yend = PL_Dpr_mean),
               color = "black",
               linewidth = 0.75) +
  ggtitle(Dprvals) +
  ylab("D'") +
  theme(legend.position="none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"))

# r^2
rsqplot <- ggplot(HZLDscores, 
                  aes(x=Same_Chr, 
                      y=r_sq, 
                      color=Same_Chr)) +
  geom_jitter(position=position_jitter(0.2),
              cex=3,
              alpha=0.75) +
  scale_x_discrete(labels = c("Physically Unlinked", 
                              "Physically Linked")) +
  geom_segment(aes(x = 0.75,
                   y = PUL_rsq_mean,
                   xend = 1.25,
                   yend = PUL_rsq_mean),
               color = "black",
               linewidth = 0.75) +
  geom_segment(aes(x = 1.75,
                   y = PL_rsq_mean,
                   xend = 2.25,
                   yend = PL_rsq_mean),
               color = "black",
               linewidth = 0.75) +
  geom_text(aes(x = 2,
                y = max(c0pHZLD$`R^2`, na.rm = T),
                label = "c0p356-c0p628"),
            size = 3,
            color = "black",
            nudge_y = -0.03) +
  ggtitle(rsqvals) +
  ylab("r^2") +
  theme(legend.position="none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, 
                                   linetype = "solid", 
                                   colour = "black"))

Allplots <- ggarrange(Dplot, DPplot, rsqplot,
                      ncol = 3, nrow = 1)
png(filename = "LDplots15jun23FINAL.png",
    width = 10,
    height = 6,
    units = "in",
    res = 1000)
Allplots
dev.off()


nogen <- c()
for (i in 1:length(colnames(c0pHZGens))) {
  nogen[i] <- unlist(strsplit(colnames(c0pHZGens)[i],
                              split = "_"))[1]
}
c0pHZGens2 <- c0pHZGens
colnames(c0pHZGens2) <- nogen
c0pHZLD2 <- LD(c0pHZGens2)

library(reshape2)
c0pHZLD_D <- melt(c0pHZLD2$D,
                  na.rm = T)
c0pHZLD_Dpr <- melt(c0pHZLD2$`D'`,
                  na.rm = T)
c0pHZLD_rsq <- melt(c0pHZLD2$`R^2`,
                  na.rm = T)




# Heatmaps
library(ggplot2)
# D 
dheat <- ggplot(data=c0pHZLD_D, 
       aes(Var2, Var1, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="darkred", high="darkblue", mid="white", midpoint = 0.1, space="Lab", name="Pairwise D") +
  xlab("") +
  ylab("") +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, 
                                 size=10, hjust=1)) +
  coord_fixed()

# D' 
dprheat <- ggplot(data=c0pHZLD_Dpr, 
                aes(Var2, Var1, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="darkred", high="darkblue", mid="white", midpoint = 0.5, space="Lab", name="Pairwise D'") +
  xlab("") +
  ylab("") +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, 
                                 size=10, hjust=1)) +
  coord_fixed()

# r^2 
rsqheat <- ggplot(data=c0pHZLD_rsq, 
                aes(Var2, Var1, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="darkred", high="darkblue", mid="white", midpoint = 0.2, space="Lab", name="Pairwise r^2") +
  xlab("") +
  ylab("") +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, 
                                 size=10, hjust=1)) +
  coord_fixed()
heatplots <- ggarrange(dheat, dprheat, rsqheat,
                      ncol = 1, nrow = 3)
setwd(genWDlaptop)
png(filename = "LDheatmaps.png",
    width = 6,
    height = 10,
    units = "in",
    res = 1000)
heatplots
dev.off()
