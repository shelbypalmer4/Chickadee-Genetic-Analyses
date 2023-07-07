#### dealing with STRUCTURE output ####

strWDmac <- "/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses/HZCH_01May2023"
strWDlaptop <- "C:/Users/Shelby Palmer/Desktop/CHICKADEES/Chickadee-Genetic-Analyses/HZCH_01May2023"

genWDmac <- "/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses"
genWDlaptop <- "C:/Users/Shelby Palmer/Desktop/CHICKADEES/Chickadee-Genetic-Analyses"

setwd(strWDmac)
setwd(strWDlaptop)
data <- read.table("HZCH_STRUCTURE_LM.output")

setwd(genWDlaptop)
setwd(genWDmac)
metadata <- read.table("HZCH_StructureInput_Final.txt")

ind <- unique(metadata[,1])

pops <- metadata$V2
pops <- pops[seq(1,length(pops),2)]

data <- data.frame(ind_ID = ind,
                   pop_ID = pops,
                   prob_CA = data$V6,
                   prob_BC = data$V7)

# this is fine for starters, but I need to be able to manipulate the axis labels and add a legend
# library(conStruct)
data <- data[order(data$pop_ID),]
# data1 <- as.matrix(data)
# make.structure.plot(admix.proportions = data1[,3:4],
#                     mar = c(4,4,2,2),
#                     sample.names = data1[,1],
#                     layer.colors = c("red4", "cyan3"))

# this looks pretty good
library(ggplot2)
ggplot(data) +
  geom_col(aes(x=ind_ID, y=1, fill="red")) +
  geom_col(aes(x=ind_ID, y=prob_CA, fill="blue")) +
  xlab("Individual ID") +
  ylab("Probability of assignment (Q)") +
  scale_fill_discrete(name = "Population", labels = c("1","2")) +
  scale_x_discrete(limits = data$ind_ID) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())

# can fine-tune the order of the bars further by adding coordinates to the dataframe with the Q-values
setwd(genWDmac)
write.csv(data, "STRUCTURE_data_with_localities.csv")
data_loc <- read.csv("STRUCTURE_data_with_localities.csv")

data_loc <- data_loc[order(data_loc$latitude),]
# this isn't as simple as I'd anticipted, since the HZ really runs southwest to northeast over here; the BCCH Bates CO location is actually barely south of SPPUA

#
#
#
# 26-27 May 2023
#### Generating categories for ancestry ####
# roughly following McQuillan et al. 2018 Van Huynh and Rice 2020

setwd(genWDmac)
gen <- read.csv("STRUCTURE_data_with_localities.csv")
singerIDs <- c("HZ3", "HZ16", "HZ20", "HZ22", "HZ24", "HZ25", "HZ27", "HZ28", "HZ30", "HZ36")

# look at values for fully-recorded singers
hist(gen$prob_CA[which(gen$ind_ID %in% singerIDs)])
gen[which(gen$ind_ID %in% singerIDs),]

# add locality-based species initials to the dataframe
GetSp <- function(x) {
  paste(unlist(strsplit(x, split = ""))[1:2], collapse = "")
}
GetSp(gen$ind_ID[34]) # works
Sp <- unlist(lapply(gen$ind_ID, GetSp))
gen <- cbind(gen, Sp)

# generate 90% CI for known CACH and BCCH
CAgens <- gen$prob_CA[which(gen$Sp == "CA")]
BCgens <- gen$prob_CA[which(gen$Sp == "BC")]

CI <- function(x, alpha) {
  t_score = qt(p=alpha/2, 
               df=length(x)-1,
               lower.tail=F)
  margin_error <- t_score * sd(x)/sqrt(length(x))
  CIL <- mean(x) - margin_error
  CIU <- mean(x) + margin_error
  print(c(CIL, CIU))
}
CA_CI <- CI(CAgens, 0.10)
# [1] 0.924495 0.989705

BC_CI <- CI(BCgens, 0.10)
# [1] 0.01187402 0.03339871

# add CI-based species initials to the frame
gen$Sp_assigned <- c(rep(NA))
gen$Sp_assigned[which(gen$prob_CA>CA_CI[1])] <- "CA"
gen$Sp_assigned[which(gen$prob_CA<BC_CI[2])] <- "BC"
gen$Sp_assigned[which(gen$prob_CA<CA_CI[1] & gen$prob_CA>BC_CI[2])] <- "HY"

# cleaning up
colnames(gen)[8] <- "Sp_by_locality"
gen <- replace(gen$Sp_by_locality, 
               which(gen$Sp_by_locality=="HZ"), 
               "HY")

# write a new csv with info needed for ancestry categorization
write.csv(gen, "STRUCTURE_data_with_localities.csv")

citation(package = "genetics")

# To cite package ‘genetics’ in publications use:
#   
#   Warnes G, Gorjanc wcfG, Leisch F, Man. M (2021). _genetics:
#   Population Genetics_. R package version 1.3.8.1.3,
# <https://CRAN.R-project.org/package=genetics>.

# summary stats for Sparrowfoot Q scores
mean(data$prob_CA[which(data$pop_ID==7)])
sd(data$prob_CA[which(data$pop_ID==7)])

# summary stats for Clinton Q scores
mean(data$prob_CA[which(data$pop_ID==6)])
sd(data$prob_CA[which(data$pop_ID==6)])

# summary stats for BCCH Q scores
mean(data$prob_CA[which(data$pop_ID %in% 1:5)])
sd(data$prob_CA[which(data$pop_ID %in% 1:5)])

# summary stats for CACH Q scores
mean(data$prob_CA[which(data$pop_ID %in% 8:10)])
sd(data$prob_CA[which(data$pop_ID %in% 8:10)])
