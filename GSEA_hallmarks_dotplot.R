setwd("~/Documents/makeBubblePlotFromGSEAstats/")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(stringr)
library(matrixStats)


dataDir = "/Users/kaalexa/gsea_home/output/jun30"

# prefix of the folders
prefix = "HALLMARK"

# cancer list
#CANCERS = c("KIRC", "SKCM", "BRCA")
#CANCERS = c("KIRC", "SKCM", "BRCA",  "LIHC", "PRAD")
CANCERS = c("KIRC", "SKCM", "BRCA", "BLCA", "CESC", "COAD", "ESCA", "HNSC", "KIRP", "LIHC", "PAAD", "PCPG", "PRAD", "SARC", "STAD", "THCA", "THYM", "UCEC")

# first cancer is KIRC
cancer = CANCERS[1]
dir = paste(dataDir, "/", prefix, "_", cancer, ".GseaPreranked.", sep = "")
dir = list.dirs(path = dataDir, full.names = TRUE)[grep(dir, list.dirs(path = dataDir, full.names = TRUE))[1]]
negFile = list.files(path = dir, pattern = "gsea_report_for_na_neg_[0-9]+.tsv", full.names = TRUE)
posFile = list.files(path = dir, pattern = "gsea_report_for_na_pos_[0-9]+.tsv", full.names = TRUE)
neg <- read.table(negFile, sep = "\t", header=T)
pos <- read.table(posFile, sep = "\t", header=T)
MERGED <- rbind(neg,pos)
MERGED$CANCER <- "1KIRC"

# Assemble data
for (cancer in CANCERS[2:length(CANCERS)]){
  dir = paste(dataDir, "/", prefix, "_", cancer, ".GseaPreranked.", sep = "")
  dir = list.dirs(path = dataDir, full.names = TRUE)[grep(dir, list.dirs(path = dataDir, full.names = TRUE))[1]]
  negFile = list.files(path = dir, pattern = "gsea_report_for_na_neg_[0-9]+.tsv", full.names = TRUE)
  posFile = list.files(path = dir, pattern = "gsea_report_for_na_pos_[0-9]+.tsv", full.names = TRUE)
  neg <- read.table(negFile, sep = "\t", header=T)
  pos <- read.table(posFile, sep = "\t", header=T)  
  CURRENT <- rbind(neg, pos)
  CURRENT$CANCER <- cancer
  MERGED <- rbind(MERGED, CURRENT)
}

MERGED$negLOGq <- -log10(MERGED$FDR.q.val)
MERGED$negLOGq[MERGED$negLOGq == "Inf"] <- 5

# Take terms significant in any cancer type
pTable <- MERGED %>% group_by(CANCER,NAME) %>% summarise(FDR.q.val) %>% pivot_wider(names_from = CANCER, values_from = FDR.q.val)
#toKeep <- pTable$NAME[rowMins(as.matrix(pTable[,2:length(pTable)])) < 0.005]

#Terms significant in KIRC
toKeep <- MERGED$NAME[MERGED$CANCER == "1KIRC" & MERGED$FDR.q.val < 0.05]


SUBSET <- MERGED[MERGED$NAME %in% toKeep,]
  


p = SUBSET%>%
  ggplot(aes(x=CANCER,y=reorder(NAME, -NES))) +
  geom_point(aes(size= negLOGq,color= NES) ) +
  #next line add border
  geom_point(aes(size= negLOGq),shape = 1,colour = "black") +
  scale_color_gradient2(low = "#27495C", mid = "white", high = "#C66B3C", midpoint = 0)+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 10),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 7), breaks = seq(0, 5, by = 1), limits = c(0,5)) +
  #reverse Y
  scale_y_discrete(limits = rev)

filename = paste(prefix, "_dotplot.pdf", sep = "")
pdf(filename, width = 8.5, height = 5, onefile=FALSE)
print(p)
dev.off()