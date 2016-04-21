# Author: Ioannis Moustakas, i.moustakas@uva.nl (Based on Martijs Jonker DEGs Script, m.j.jonker@uva.nl )
# Title: Perform a DEG analysis with DESeq2 R package
# Save the comparison pairs into a table and feed it as input. This is handy when a large number of comparisons is required.
# There is an annotation function that is serving to anontate the gene names with a short description of the product function  

library(DESeq2)

setwd("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/")
options(stringsAsFactors = T)

pathDesign=c("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/")
# load the design file 
design <- read.delim(paste0(pathDesign,"designDEGAnalysis.csv"),stringsAsFactors=FALSE)
designFile <- design
row.names(designFile) <- design$Sample 

pathData="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/mapping/"
# load the gene count table
countTable <- read.delim(paste0(pathData, "allSamplesCountTable.txt"),stringsAsFactors=FALSE)
nrows <- nrow(countTable) 
ncolns <- ncol(countTable)
countTable <- countTable[-c((nrows-4):nrows), ] 
geneIDs <- as.vector(countTable[,1])
countData <- data.matrix(countTable[,2:ncolns])
row.names(countData) <- geneIDs

# Build the countTable and designFile for Block 1: Samples S01-S05
blockOneCountTable <- countData[ ,1:5]
blockOneDesignFile <- designFile[1:5,]

############ Try Out DESeq1 ############
library(DESeq)

cds = newCountDataSet(blockOneCountTable, blockOneDesignFile$SampleID)
cds <- estimateSizeFactors(cds)
cds2 <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
res2 <- nbinomTest(cds2, "eTD", "e3p")
plotDispEsts(cds2)




#################### $$$$$$$$$$$$$$$$$ #################

dds <- DESeqDataSetFromMatrix(countData = blockOneCountTable, colData = blockOneDesignFile, design = ~ SampleID )
 



res <- results(dds, contrast = c("SampleID", "eTD", "e3n"))
plotMA(res, main=fileName, ylim = c(-5, 5))
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

# Build the DESeq data set object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = designFile, design = ~ SampleID )
dds <- DESeq(dds)

# Load the table with the comparissons to make
allComparissons <- read.delim(paste0(pathDesign, "comparissonsMatrix.csv"), stringsAsFactors = F)

# create data.frame to store some sample statistics
statistics <- data.frame(comparisson=c(1), pairContrast=c(1), minPadj=c(1), numGenesPadjBelThres=c(1))

# Where the gene tables are saved
tablePath = paste0(pathDesign, "tables/")
dir.create(tablePath)

# Where the plots are saved
plotPath = paste0(pathDesign, "plots/")
dir.create(plotPath)

# Do the comparissons as they are defined in the comparisson table
for (i in 1:nrow(allComparissons)) {
  x <- allComparissons[i, ]
  fileName = paste(x[,2], x[,3], sep = ".VS.")
  
  # create Dir where the results are saved
  dir.create( paste0(plotPath, x[,1]), showWarnings = F)
    
  # define the pair of samples to contrast
  res <- results(dds, contrast = c("SampleID", x[,2], x[,3]))
  
  completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
  png(completePlotPath, width=635,height=460)
  plotMA(res, main=fileName, ylim = c(-5, 5))
  graphics.off()
  
  minPadj <- min(res$padj[!is.na(res$padj)])
  geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
  comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
  statistics <- rbind(statistics, comparissonStats)
  
  # order the table on the padj values, so the most significant genes are on top
  orderedTable <- res[order(res$padj), ]
  
  dir.create( paste0(tablePath, x[,1]), showWarnings = F)
  completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
  write.table(orderedTable, completeTablePath, sep="\t")
  
}

#################### $$$$$$$$$$$$$$$$$$$$$$$$ ########################
#### Need to use another column of the design table do the last 3 comparissons of set1
dds <- DESeqDataSetFromMatrix(countData = countData, colData = designFile, design = ~ Sample )
dds <- DESeq(dds)

########
# Do the S07 VS S17
fileName= "f3-24.VS.f304-24"
# define the pair of samples to contrast
res <- results(dds, contrast = c("Sample", "S07", "S17"))

completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
png(completePlotPath, width=635,height=460)
plotMA(res, main=fileName, ylim = c(-5, 5))
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
statistics <- rbind(statistics, comparissonStats)

# order the table on the padj values, so the most significant genes are on top
orderedTable <- res[order(res$padj), ]

dir.create( paste0(tablePath, x[,1]), showWarnings = F)
completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
write.table(orderedTable, completeTablePath, sep="\t")

########
# Do the S07 VS S18
fileName= "f3-24.VS.f314-24"
# define the pair of samples to contrast
res <- results(dds, contrast = c("Sample", "S07", "S18"))

completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
png(completePlotPath, width=635,height=460)
plotMA(res, main=fileName, ylim = c(-5, 5))
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
statistics <- rbind(statistics, comparissonStats)

# order the table on the padj values, so the most significant genes are on top
orderedTable <- res[order(res$padj), ]

dir.create( paste0(tablePath, x[,1]), showWarnings = F)
completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
write.table(orderedTable, completeTablePath, sep="\t")

########
# Do the S07 VS S19
fileName= "f3-24.VS.f315-24"
# define the pair of samples to contrast
res <- results(dds, contrast = c("Sample", "S07", "S19"))

completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
png(completePlotPath, width=635,height=460)
plotMA(res, main=fileName, ylim = c(-5, 5))
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
statistics <- rbind(statistics, comparissonStats)

# order the table on the padj values, so the most significant genes are on top
orderedTable <- res[order(res$padj), ]

dir.create( paste0(tablePath, x[,1]), showWarnings = F)
completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
write.table(orderedTable, completeTablePath, sep="\t")

# Write statistics TAble to a file
completeTablePath = paste0(tablePath, "cumulativeStatistics.txt")
write.table(statistics, completeTablePath, sep="\t", row.names = F)
