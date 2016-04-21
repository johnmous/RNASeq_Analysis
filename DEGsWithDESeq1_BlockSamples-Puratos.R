# Author: Ioannis Moustakas, i.moustakas@uva.nl (Based on Martijs Jonker DEGs Script, m.j.jonker@uva.nl )
# Title: Perform a DEG analysis with DESeq R package
# Save the comparison pairs into a table and feed it as input. This is handy when a large number of comparisons is required.
# There is an annotation function that is serving to anontate the gene names with a short description of the product function  

library(DESeq)
library(limma)

compare <- function(block){
  for (i in 1:nrow(block)) {
    # i=1
    x <- block[i, ]
    fileName = paste(x[,2], x[,3], sep = ".VS.")
    res <- nbinomTest(cds2, x[,2], x[,3])
    
    # annotate BSU number 
    BSUs <- res$id
    productDescr <- vector()
    for (BSU in BSUs) {
      productDescr <- c(productDescr, annotTable[annotTable$BSU==BSU, ]$productVector)
    }
    
    dir.create(paste0(plotPath, x[,1]), showWarnings = F)
    
    # plot and save the MA plot
#     completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
#     png(completePlotPath, width=1280,height=1024)
#     plotMA(res)
#     graphics.off()
    
    # Calculate some statistics on about the sample and save them
    minPadj <- min(res$padj[!is.na(res$padj)])
    geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
    comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
    statistics <- rbind(statistics, comparissonStats)
    
    ## plot Padj Histogram
#     completePlotPath <- paste0(plotPath, x[,1],"/histPadj_", fileName, ".png")
#     png(completePlotPath, width=1280,height=1024)
#     hist(res$pval,breaks=30)
#     graphics.off()
    
    # Add an annotation (description) column  
    res$productDescr <- productDescr
    # Re-Order the data frame to bring the description next to geneID
    res <- res[, c(1,9,2:8)] 
    dir.create(paste0(tablePath, x[,1]), showWarnings = F)
    ## create Dir
    completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
    write.table(res, completeTablePath, sep="\t", row.names = F)
    
  }
}

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

# create data.frame to store some sample statistics
statistics <- data.frame(comparisson=c(1), pairContrast=c(1), minPadj=c(1), numGenesPadjBelThres=c(1))

# load the annotation table (build from the gff file)
annotTable <- read.delim(paste0(pathDesign, "AnnotationAppendix.txt"),stringsAsFactors=FALSE)

outputPath="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/DESeqOne/"
# Where the gene tables are saved
tablePath = paste0(outputPath, "tables/")
dir.create(tablePath, showWarnings = F)

# Where the plots are saved
plotPath = paste0(outputPath, "plots/")
dir.create(plotPath, showWarnings = F)

############ $$$$$$$$$$ ###########
#Build the countTable and designFile for Block 1: Samples S01-S05
blockOneCountTable <- countData[ ,1:5]
blockOneDesignFile <- designFile[1:5,]

# DESeq1
cds = newCountDataSet(blockOneCountTable, blockOneDesignFile$SampleID)
cds <- estimateSizeFactors(cds)
cds2 <- estimateDispersions(cds, method="blind", sharingMode="fit-only")



completePlotPath= paste0(plotPath, "/DispEstim_BlockOne.png")
png(completePlotPath, width=635,height=460)
plotDispEsts(cds2)
graphics.off()

blockOneComparissons <- data.frame(Set=c("Block1", "Block1", "Block1", "Block1"),
                              GroupOne=c("eTD", "eTD", "eTD", "eTD"),
                              GroupTwo=c("e3p", "e3n", "e4p", "e4n"), stringsAsFactors = F)

compare(blockOneComparissons)

############ $$$$$$$$$$ ###########
#Build the countTable and designFile for Block 2: Samples S06-S09
blockTwoCountTable <- countData[ ,6:9]
blockTwoDesignFile <- designFile[6:9, ]

# DESeq1
cds = newCountDataSet(blockTwoCountTable, blockTwoDesignFile$SampleID)
cds <- estimateSizeFactors(cds)
cds2 <- estimateDispersions(cds, method="blind", sharingMode="fit-only")



completePlotPath= paste0(plotPath, "/DispEstim_BlockTwo.png")
png(completePlotPath, width=635,height=460)
plotDispEsts(cds2)
graphics.off()

blockTwoComparissons <- data.frame(Set=c("Block2", "Block2", "Block2"),
                                   GroupOne=c("fTD-24", "fTD-24", "fTD-24"),
                                   GroupTwo=c("f3-24", "f4A-24", "f4B-24"), stringsAsFactors = F)

compare(blockTwoComparissons)

############ $$$$$$$$$$ ###########
#Build the countTable and designFile for Block 3: Samples S07 and S17-19 
blockThreeCountTable <- countData[ ,c(7, 17:19)]
blockThreeDesignFile <- designFile[c(7, 17:19), ]

# DESeq1
cds = newCountDataSet(blockThreeCountTable, blockThreeDesignFile$Sample)
cds <- estimateSizeFactors(cds)
cds2 <- estimateDispersions(cds, method="blind", sharingMode="fit-only")

completePlotPath= paste0(plotPath, "/DispEstim_BlockThree.png")
png(completePlotPath, width=635,height=460)
plotDispEsts(cds2)
graphics.off()

blockThreeComparissons <- data.frame(Set=c("Block3", "Block3", "Block3"),
                                   GroupOne=c("S07", "S07", "S07"),
                                   GroupTwo=c("S17", "S18", "S19"), stringsAsFactors = F)

compare(blockThreeComparissons)













# ############ $$$$$$$$$$ ###########
# #Build the countTable and designFile for Set 2: Samples S14-25 
# blockFourCountTable <- countData[ ,14:25]
# blockFourDesignFile <- designFile[14:25, ]
# #blockFourDesignFile <- blockFourDesignFile[,c(4,5)]
# 
# # DESeq1
# cds = newCountDataSet(blockFourCountTable, blockFourDesignFile)
# cds <- estimateSizeFactors(cds)
# #cds2 <- estimateDispersions(cds, method="blind",sharingMode="fit-only" )
# cds2 <- estimateDispersions(cds, method="pooled-CR", modelFormula = count ~ SampleID + Culture )
# 
# plotDispEsts(cds2)
# 
# # Fit a generalized linear model (GLM) for each gene. 
# # SampleID is the time point: 3h (exp), 24h, 48h, 65h. 
# # Culture is the fermentor ID.
# fit1 <- fitNbinomGLMs(cds2, count ~ SampleID + Culture)
# fit0 <- fitNbinomGLMs(cds2, count ~ Culture )
# 
# # Perform chi-squared tests comparing two sets of GLM fits
# pvalsGLM <- nbinomGLMTest(fit1, fit0)
# hist(pvalsGLM, breaks = 100)
# padjGLM <- p.adjust(pvalsGLM, method="BH")
# 
# # Number of DEGs, when taking into account only the effect of SampleId (time points) removing the effect of Culture (fermentor)
# nOfDEGs <- length(which(padjGLM<0.05))
# 
# # save the table with all differentially expressed genes
# DEGsBlockFour <- blockFourCountTable[which(padjGLM<0.05), ]
# 
# # A vector with all DEGs names (BSUs)
# BSUsDEGsBlockFour <- row.names(blockFourCountTable[which(padjGLM<0.05), ])
# 
# 
# ######################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ##############################
# # Pathway analysis. Use the code from Paul. 
# ##### Kegg pathways
# # R-3.0.1
# # Use hypergeometric test 
# 
# biocLite("KEGGREST")
# library("KEGGREST")
# library("multtest")
# 
# listDatabases()
# ls("package:KEGGREST") 
# 
# 
# # bsu genes
# bsu.genes <- keggList("bsu")
# bsu.genes <- names(bsu.genes)
# length(unique(bsu.genes)) # [1] 4421
# 
# # bsu pathways
# bsu.pw <- keggList("pathway", "bsu")
# length(bsu.pw) # [1] 114
# 
# # Get genes per pathway
# genes.pw <- list()
# for(i in 1:length(bsu.pw)){
#   temp <- keggGet(names(bsu.pw[i]))
#   genes.temp <-  temp[[1]]$GENE
#   if(length(genes.temp) > 0){
#     genes.temp <- genes.temp[seq(1, length(genes.temp), by=2)]
#     genes.pw[[names(bsu.pw[i])]] <- paste("bsu:", genes.temp, sep="")
#     #genes.pw[[names(bsu.pw[i])]] <- genes.temp
#   }
# }
# length(genes.pw) # [1] 102
# 
# DEGs <- row.names(DEGsBlockFour)
# DEG.names <- paste0("bsu:", DEGs) 
# 
# # Over representation
# hyper.test <- function(gene.set, selected, gene.list)
#   # gene.set = pathway genes
#   # selected = DEGs
#   # gene.list = Background
# {
#   numWdrawn <- length(intersect(gene.set, selected))
#   numDrawn <- length(selected)
#   numW <- length(gene.set)
#   numB <- length(gene.list) - numW
#   
#   #over represented
#   pvalue <- phyper(numWdrawn - 1, numW, numB, numDrawn, lower.tail=FALSE) #P(X > x)
#   
#   #under represented
#   #pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=TRUE) #P(X <= x)
#   
#   x <- max(0, numDrawn-numB):min(numW, numDrawn)
#   y <- dhyper(x, numW, numB, numDrawn)
#   names(y) <- x
#   #return(list(observed=numWdrawn, estimated=y, pvalue=pvalue))
#   return(pvalue)
# } 
# 
# 
# pw <- genes.pw
# Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=bsu.genes))
# pvals <- as.vector(Overrep)
# T.res <- mt.rawp2adjp(pvals, "BH")
# adjPvals <- T.res$adjp[order(T.res$index),]
# rownames(adjPvals) <- names(Overrep)
# Overrep <- names(Overrep[which(Overrep < 0.1)])
# Overrep_adj <- adjPvals[as.vector(which(adjPvals[,2] < 0.1)),] 
# 
# 
# ##### Kegg pathways 
# # Use geneSetTest (based on Wilcoxon signed-rank test)
# 
# bsu.pw <- keggList("pathway", "bsu")
# length(bsu.pw) # [1] 114
# 
# # Get genes per pathway
# genes.pw <- list()
# for(i in 1:length(bsu.pw)){
#   temp <- keggGet(names(bsu.pw[i]))
#   genes.temp <-  temp[[1]]$GENE
#   if(length(genes.temp) > 0){
#     genes.temp <- genes.temp[seq(1, length(genes.temp), by=2)]
#     #genes.pw[[names(bsu.pw[i])]] <- paste("bsu:", genes.temp, sep="")
#     genes.pw[[names(bsu.pw[i])]] <- genes.temp
#   }
# }
# length(genes.pw) # [1] 102
#  
# FC <- DEGsBlockFour
# 
# # prepare geneSetTest
# gst.result <- list()
# for(i in 1:length(genes.pw)){
#   m <- match(genes.pw[[i]], rownames(FC))
#   rm <- which(is.na(m) == TRUE)
#   if(length(rm) > 0){
#     m <- m[-rm]
#   }
#   genes.get <- as.vector(rownames(FC)[m], mode="character")
#   gst.temp <- c()
#   for(j in 1:ncol(FC)){
#     gst.temp <- c(gst.temp, geneSetTest(index=m, statistics=FC[,j] ))
#   }
#   gst.result[[names(genes.pw)[i]]] <- gst.temp
# } 
# 
# out <- NULL
# for(i in 1:length(gst.result)){
#   Kegg <- names(bsu.pw[match(names(gst.result)[i], names(bsu.pw))])
#   pwn <- bsu.pw[[match(names(gst.result)[i], names(bsu.pw))]]
#   out <- rbind(out, c(Kegg, pwn, gst.result[[i]]))
# }
# 
# colnames(out) <- c("KeggID","Pathway", designFile$SampleID[14:25] )
# head(out)
# write.table(out, paste0(pathDesign, "KEGG/Pathway_geneSetTest.txt"), quote=F, sep="\t", col.names=T, row.names=F)             
# 
# 
# 
# ##### $$$$$ #####
# # Compile a list with all BSUs per GO term
# 
# BsubtilisGOA = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/6.B_subtilis_168.goa.txt"
# GOATable <- read.delim(BsubtilisGOA, header = F, comment.char = "#")
# 
# # Build BSUs vector 
# # First Separate on pipe ("|"), grep only the elements containing BSUs and separate again on "/"
# BSUs <- sapply(strsplit(as.character(GOATable$V11), "|", fixed=T), function(x) {
#         index <- grep("BSU", x, fixed = T)
#         #x[index]
#         unlist(strsplit(x[index], "/", fixed=T))
#       })
# 
# # sapply(BSUs, function(x){
# #   length(x)
# # })
# 
# # extract all go terms present in bacillus annotation table and asign them to the BSUs list. 
# # This way a GO => BSU mapping is build out of the annotation table
# GOTerms <- GOATable$V5
# names(BSUs) <- GOTerms
# 
# # aggregate all BSU that map to a GO term together so we build a GO => [BSUs] 
# allGoTerms <- unique(GOTerms)
# GOToBSU <- NULL
# i = 1
# for (term in allGoTerms){
#   BSUVector <- unique(unlist(BSUs[ names(BSUs) == term ]))
#   GOToBSU[[i]] <- BSUVector
#   i = i + 1
# }
# names(GOToBSU) <- allGoTerms
# 
# ####### GO enrichement analysis 
# FC <- DEGsBlockFour
# genes.pw <- GOToBSU
# 
# # prepare geneSetTest
# gst.result <- list()
# for(i in 1:length(genes.pw)){
#   m <- match(genes.pw[[i]], rownames(FC))
#   rm <- which(is.na(m) == TRUE)
#   if(length(rm) > 0){
#     m <- m[-rm]
#   }
#   genes.get <- as.vector(rownames(FC)[m], mode="character")
#   gst.temp <- c()
#   for(j in 1:ncol(FC)){
#     gst.temp <- c(gst.temp, geneSetTest(index=m, statistics=FC[,j] ))
#   }
#   gst.result[[names(genes.pw)[i]]] <- gst.temp
# }
# 
# out <- NULL
# for(i in 1:length(gst.result)){
#   Kegg <- names(bsu.pw[match(names(gst.result)[i], names(bsu.pw))])
#   pwn <- bsu.pw[[match(names(gst.result)[i], names(bsu.pw))]]
#   out <- rbind(out, c(Kegg, pwn, gst.result[[i]]))
# }
# 
# 
# length(unique(GOToBSU[[1]]))
# 
# 
# ddply()
# 
# length(unique(BSUs))
# length(BSUs)
# length(unlist(BSUs))
# 
# 
# length(unique(GOTerms))
# length(GOTerms)
# 
# lapply(strsplit(as.character(GOATable$V11), "|", fixed=T), function(x) {
#   index <- grep("BSU", x)
# })
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #################### $$$$$$$$$$$$$$$$$ #################
# 
# dds <- DESeqDataSetFromMatrix(countData = blockOneCountTable, colData = blockOneDesignFile, design = ~ SampleID )
#  
# 
# 
# 
# res <- results(dds, contrast = c("SampleID", "eTD", "e3n"))
# plotMA(res, main=fileName, ylim = c(-5, 5))
# minPadj <- min(res$padj[!is.na(res$padj)])
# geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
# 
# # Build the DESeq data set object
# dds <- DESeqDataSetFromMatrix(countData = countData, colData = designFile, design = ~ SampleID )
# dds <- DESeq(dds)
# 
# # Load the table with the comparissons to make
# allComparissons <- read.delim(paste0(pathDesign, "comparissonsMatrix.csv"), stringsAsFactors = F)
# 
# # create data.frame to store some sample statistics
# statistics <- data.frame(comparisson=c(1), pairContrast=c(1), minPadj=c(1), numGenesPadjBelThres=c(1))
# 
# # Where the gene tables are saved
# tablePath = paste0(pathDesign, "tables/")
# dir.create(tablePath)
# 
# # Where the plots are saved
# plotPath = paste0(pathDesign, "plots/")
# dir.create(plotPath)
# 
# # Do the comparissons as they are defined in the comparisson table
# for (i in 1:nrow(allComparissons)) {
#   x <- allComparissons[i, ]
#   fileName = paste(x[,2], x[,3], sep = ".VS.")
#   
#   # create Dir where the results are saved
#   dir.create( paste0(plotPath, x[,1]), showWarnings = F)
#     
#   # define the pair of samples to contrast
#   res <- results(dds, contrast = c("SampleID", x[,2], x[,3]))
#   
#   completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
#   png(completePlotPath, width=635,height=460)
#   plotMA(res, main=fileName, ylim = c(-5, 5))
#   graphics.off()
#   
#   minPadj <- min(res$padj[!is.na(res$padj)])
#   geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
#   comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
#   statistics <- rbind(statistics, comparissonStats)
#   
#   # order the table on the padj values, so the most significant genes are on top
#   orderedTable <- res[order(res$padj), ]
#   
#   dir.create( paste0(tablePath, x[,1]), showWarnings = F)
#   completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
#   write.table(orderedTable, completeTablePath, sep="\t")
#   
# }
# 
# #################### $$$$$$$$$$$$$$$$$$$$$$$$ ########################
# #### Need to use another column of the design table do the last 3 comparissons of set1
# dds <- DESeqDataSetFromMatrix(countData = countData, colData = designFile, design = ~ Sample )
# dds <- DESeq(dds)
# 
# ########
# # Do the S07 VS S17
# fileName= "f3-24.VS.f304-24"
# # define the pair of samples to contrast
# res <- results(dds, contrast = c("Sample", "S07", "S17"))
# 
# completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
# png(completePlotPath, width=635,height=460)
# plotMA(res, main=fileName, ylim = c(-5, 5))
# graphics.off()
# 
# minPadj <- min(res$padj[!is.na(res$padj)])
# geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
# comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
# statistics <- rbind(statistics, comparissonStats)
# 
# # order the table on the padj values, so the most significant genes are on top
# orderedTable <- res[order(res$padj), ]
# 
# dir.create( paste0(tablePath, x[,1]), showWarnings = F)
# completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
# write.table(orderedTable, completeTablePath, sep="\t")
# 
# ########
# # Do the S07 VS S18
# fileName= "f3-24.VS.f314-24"
# # define the pair of samples to contrast
# res <- results(dds, contrast = c("Sample", "S07", "S18"))
# 
# completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
# png(completePlotPath, width=635,height=460)
# plotMA(res, main=fileName, ylim = c(-5, 5))
# graphics.off()
# 
# minPadj <- min(res$padj[!is.na(res$padj)])
# geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
# comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
# statistics <- rbind(statistics, comparissonStats)
# 
# # order the table on the padj values, so the most significant genes are on top
# orderedTable <- res[order(res$padj), ]
# 
# dir.create( paste0(tablePath, x[,1]), showWarnings = F)
# completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
# write.table(orderedTable, completeTablePath, sep="\t")
# 
# ########
# # Do the S07 VS S19
# fileName= "f3-24.VS.f315-24"
# # define the pair of samples to contrast
# res <- results(dds, contrast = c("Sample", "S07", "S19"))
# 
# completePlotPath= paste0(plotPath, x[,1],"/MAPlot_", fileName, ".png")
# png(completePlotPath, width=635,height=460)
# plotMA(res, main=fileName, ylim = c(-5, 5))
# graphics.off()
# 
# minPadj <- min(res$padj[!is.na(res$padj)])
# geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
# comparissonStats <- data.frame(comparisson = x[,1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
# statistics <- rbind(statistics, comparissonStats)
# 
# # order the table on the padj values, so the most significant genes are on top
# orderedTable <- res[order(res$padj), ]
# 
# dir.create( paste0(tablePath, x[,1]), showWarnings = F)
# completeTablePath = paste0(tablePath, x[,1],"/DESeqTable_", fileName, ".txt")
# write.table(orderedTable, completeTablePath, sep="\t")
# 
# # Write statistics TAble to a file
# completeTablePath = paste0(tablePath, "cumulativeStatistics.txt")
# write.table(statistics, completeTablePath, sep="\t", row.names = F)
