#### dealing with STRUCTURE output ####

setwd("/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses/HZCH_28Apr2023")
data <- read.table("HZCH_K2.output")

setwd("/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses")
metadata <- read.table("HZCH_StructureInput_Final.txt")

ind <- unique(metadata[,1])

pops <- metadata$V2
pops <- pops[seq(1,length(pops),2)]

data <- data.frame(ind_ID = ind,
                   pop_ID = pops,
                   prob_BC = data$V6,
                   prob_CA = data$V7)

library(ggplot2)

