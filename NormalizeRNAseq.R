# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: Normalize Frank Takken RNA-seq data and check the result with the aid of ERCCs
library(ggplot2)
library(reshape2)
library(ggvis)
library(plotly)
interactive()

# load the ERCCs concetration table
concetrationTable <- read.delim("/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/ERCC/ERCC_Controls_Analysis.txt", header=T)
ERCCsConcTable <- concetrationTable[,c(2,4)]
colnames(ERCCsConcTable)[1] <- "Names"

geneAndERRCCsTable <- read.delim("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/MAPQTen/combinedERCCsGeneCountTable.txt")


# extract the ERCCs fron the geneAndERRCCsTable
ERCCsInSamples <- geneAndERRCCsTable[grepl("ERCC", geneAndERRCCsTable$Names), ]

# Normalize the geneAndERRCCsTable
# remove the bottom 5 lines that are not gene counts but the report of htseq-count
nrows <- nrow(geneAndERRCCsTable) 
countTableGenesOnly <- geneAndERRCCsTable[-c((nrows-6):nrows), ] 
# save the first columns as row name and then remove it (gene names)
row.names(countTableGenesOnly) <- countTableGenesOnly[,1]
countTableGenesOnly <- countTableGenesOnly[,-1]
# now normalize the table
# the function to normalize
sizeFactors.mad <- function (counts, locfunc = median){
  loggeomeans <- rowMeans(log(counts))
  apply(counts, 2, function(cnts) exp(locfunc((log(cnts) -
                                                 loggeomeans)[is.finite(loggeomeans)])))
}
sf <- sizeFactors.mad(countTableGenesOnly)
#divide countdata by sizefactors#
CountTable.scaled <- countTableGenesOnly
for(i in 1:ncol(CountTable.scaled)){
  CountTable.scaled[,i] <- CountTable.scaled[,i]/sf[i]
} 
write.table(CountTable.scaled,"/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/MAPQTen/NormalizedGenesERCCs.txt", sep="\t")

# extract the ERCCs fron the CountTable.scaled (Normalized)
ERCCsInNormalizedSamples <- CountTable.scaled[grepl("ERCC", row.names(CountTable.scaled)), ]
ERCCsInNormalizedSamples$Names <- row.names(ERCCsInNormalizedSamples) 
ERCCsInNormalizedSamplesConc <- merge(ERCCsInNormalizedSamples, ERCCsConcTable, by="Names")
row.names(ERCCsInNormalizedSamplesConc) <- ERCCsInNormalizedSamplesConc$Names
ERCCsInNormalizedSamplesConc[,1] <- ERCCsInNormalizedSamplesConc[,ncol(ERCCsInNormalizedSamplesConc)]
colnames(ERCCsInNormalizedSamplesConc)[1] <- "Concentration"
ERCCsInNormalizedSamplesConc <- ERCCsInNormalizedSamplesConc[,-ncol(ERCCsInNormalizedSamplesConc)]
ERCCsNormalizedCollapsed <- ddply(ERCCsInNormalizedSamplesConc, "Concentration", numcolwise(sum))
ERCCsNormalizedMelted <- melt(ERCCsNormalizedCollapsed, id="Concentration")
names(ERCCsNormalizedMelted) <- c("Concentration", "Sample", "Count")
# log trans
ERCCsNormalizedMelted$Count <- log2(ERCCsNormalizedMelted$Count+1)
ggplot(ERCCsNormalizedMelted, aes(Concentration, Count)) + ggtitle("Normalized") + scale_x_log10()+ geom_line(aes(colour = Sample))
normPlot <- qplot(Concentration, Count, data=ERCCsNormalizedMelted) + ggtitle("Normalized") + scale_x_log10()+ geom_line(aes(colour = Sample))+theme(legend.position = "right")

set_credentials_file("Ioannis.moustakas1", "ytm8z8n5em")
py <- plotly()
py$ggplotly(normPlot)

# put a column for the concetration of each of the ERCCs
ERCCsAllSamplesConc <- merge(ERCCsInSamples, ERCCsConcTable, by="Names")
# Set the name of the ERCCs as the row name
row.names(ERCCsAllSamplesConc) <- ERCCsAllSamplesConc$Names
ERCCsAllSamplesConc[,1] <- ERCCsAllSamplesConc[,ncol(ERCCsAllSamplesConc)]
colnames(ERCCsAllSamplesConc)[1] <- "Concentration"
ERCCsAllSamplesConc <- ERCCsAllSamplesConc[,-ncol(ERCCsAllSamplesConc)]
ERCCsCollapsed <- ddply(ERCCsAllSamplesConc, "Concentration", numcolwise(sum))
ERCCsMelted <- melt(ERCCsCollapsed, id="Concentration")
names(ERCCsMelted) <- c("Concentration", "Sample", "Count")
ERCCsMelted$Count <- log2(ERCCsMelted$Count+1)
ggplot(ERCCsMelted, aes(Concentration, Count)) +  ggtitle("Original") + scale_x_log10()+ geom_line(aes(colour = Sample))



########## Normalize on ERCCs only ##########

# save the first columns as row name and then remove it (gene names)
row.names(ERCCsInSamples) <- ERCCsInSamples[,1]
ERCCsInSamples <- ERCCsInSamples[,-1]

sf <- sizeFactors.mad(ERCCsInSamples)
#divide countdata by sizefactors#
CountTable.scaled <- ERCCsInSamples
for(i in 1:ncol(CountTable.scaled)){
  CountTable.scaled[,i] <- CountTable.scaled[,i]/sf[i]
} 
write.table(CountTable.scaled,"/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/MAPQTen/NormalizedOnERCCs.txt", sep="\t")


ERCCsInNormalizedSamples <- CountTable.scaled[grepl("ERCC", row.names(CountTable.scaled)), ]
ERCCsInNormalizedSamples$Names <- row.names(CountTable.scaled) 
ERCCsInNormalizedSamplesConc <- merge(ERCCsInNormalizedSamples, ERCCsConcTable, by="Names")
row.names(ERCCsInNormalizedSamplesConc) <- ERCCsInNormalizedSamplesConc$Names
ERCCsInNormalizedSamplesConc[,1] <- ERCCsInNormalizedSamplesConc[,ncol(ERCCsInNormalizedSamplesConc)]
colnames(ERCCsInNormalizedSamplesConc)[1] <- "Concentration"
ERCCsInNormalizedSamplesConc <- ERCCsInNormalizedSamplesConc[,-ncol(ERCCsInNormalizedSamplesConc)]
ERCCsNormalizedCollapsed <- ddply(ERCCsInNormalizedSamplesConc, "Concentration", numcolwise(sum))
ERCCsNormalizedMelted <- melt(ERCCsNormalizedCollapsed, id="Concentration")
names(ERCCsNormalizedMelted) <- c("Concentration", "Sample", "Count")
# log trans
ERCCsNormalizedMelted$Count <- log2(ERCCsNormalizedMelted$Count+1)
ggplot(ERCCsNormalizedMelted, aes(Concentration, Count)) + ggtitle("NormalizedOnERCCs") + scale_x_log10()+ geom_line(aes(colour = Sample))










all_values <- function(x) {
  if(is.null(x)) return(NULL)
  paste0(names(x), ": ", format(x), collapse = "<br />")
}

ERCCsNormalizedMelted %>% ggvis(~Concentration, ~Count, size.hover := 200) %>% 
  scale_numeric("x", trans="log", expand=0) %>% 
  layer_points(fill=~factor(Sample)) 

ERCCsMelted %>% ggvis(~Concentration, ~Count) %>% 
  scale_numeric("x", trans="log", expand=0) %>% 
  layer_points(fill=~factor(Sample)) %>% 
  layer_smooth(method = "lm")




# melt the ERCCs
ERCCsAllSamplesConcMelted <- melt(ERCCsAllSamplesConc[,1:5], id="Concentration")
names(ERCCsAllSamplesConcMelted) <- c("Concentration", "Sample", "Count")
# log trans
ERCCsAllSamplesConcMelted$Count <- log2(ERCCsAllSamplesConcMelted$Count+1)


all_values <- function(x) {
  if(is.null(x)) return(NULL)
  paste0(names(x), ": ", format(x), collapse = "<br />")
}

ERCCsAllSamplesConcMelted %>% ggvis(~Concentration, ~Count, size.hover := 200) %>% 
  scale_numeric("x", trans="log", expand=0) %>% layer_points(fill=~factor(Sample)) 


# Data.frame with S01 original and normalized
t <- data.frame(Concentration=ERCCsAllSamplesConc$Concentration, 
                S01=ERCCsAllSamplesConc$S01, S01Norm=ERCCsInNormalizedSamplesConc$S01)

# melt the ERCCs
ERCCsAllSamplesConcMelted <- melt(t, id="Concentration")
names(ERCCsAllSamplesConcMelted) <- c("Concentration", "Sample", "Count")