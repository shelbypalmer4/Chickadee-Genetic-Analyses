#### dealing with STRUCTURE output ####

setwd("/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses/HZCH_01May2023")
data <- read.table("HZCH_STRUCTURE_LM.output")

setwd("/Users/shelbypalmer/Documents/GitHub/Chickadee-Genetic-Analyses")
metadata <- read.table("HZCH_StructureInput_Final.txt")

ind <- unique(metadata[,1])

pops <- metadata$V2
pops <- pops[seq(1,length(pops),2)]

data <- data.frame(ind_ID = ind,
                   pop_ID = pops,
                   prob_CA = data$V6,
                   prob_BC = data$V7)

# this is fine for starters, but I need to be able to manipulate the axis labels and add a legend
library(conStruct)
data1 <- data[order(data$pop_ID),]
data1 <- as.matrix(data1)
make.structure.plot(admix.proportions = data1[,3:4],
                    mar = c(4,4,2,2),
                    sample.names = data1[,1],
                    layer.colors = c("red4", "cyan3"))

# this looks pretty good. just need to order the bars by sampling locailty
library(ggplot2)
ggplot(data) +
  geom_col(aes(x=ind_ID, y=1, fill="red")) +
  geom_col(aes(x=ind_ID, y=prob_CA, fill="blue")) +
  xlab("Individual ID") +
  ylab("Probability of assignment") +
  scale_fill_discrete(name = "Population", labels = c("1","2")) +
  scale_x_discrete(limits = data$ind_ID) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())


