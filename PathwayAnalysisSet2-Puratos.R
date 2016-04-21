# Author: Ioannis Moustakas. i.moustakas@uva.nl 
# Title: Pathway Analysis for set 2 of comparissons for Puratos RNASeq. Following directions from Martijs

library("DESeq")

# Load the relevant data
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

# load the annotation table (build from the gff file)
annotTable <- read.delim(paste0(pathDesign, "AnnotationAppendix.txt"),stringsAsFactors=FALSE)

############ $$$$$$$$$$ ###########
#Build the countTable and designFile for Set 2: Samples S14-25 
blockFourCountTable <- countData[ ,14:25]
blockFourDesignFile <- designFile[14:25, ]
#blockFourDesignFile <- blockFourDesignFile[,c(4,5)]

# DESeq1
cds = newCountDataSet(blockFourCountTable, blockFourDesignFile)
cds <- estimateSizeFactors(cds)
#cds2 <- estimateDispersions(cds, method="blind",sharingMode="fit-only" )
cds2 <- estimateDispersions(cds, method="pooled-CR", modelFormula = count ~ SampleID + Culture )

# A function to calculate the DESeq table, as calculated by nbinomTest function. 
# feed the cds object, and the sample names of the two groups that are to be contrasted
CalculateTable <- function (cds, groupA_Names, groupB_Names){
  expNames <- c(groupA_Names, groupB_Names)
  cds.part <- cds2[,expNames]
  
  fit1 <- fitNbinomGLMs(cds.part, count ~ SampleID + Culture)
  fit0 <- fitNbinomGLMs(cds.part, count ~ Culture )
  
  # Perform chi-square tests comparing two sets of GLM fits 
  pvalsGLM <- nbinomGLMTest(fit1, fit0)
  hist(pvalsGLM, breaks = 100)
  padjGLM <- p.adjust(pvalsGLM, method="BH")
  length(which(padjGLM<0.05))
  
  normCountsTable <- counts(cds.part,normalized=TRUE)
  
  ## calculate baseMean 
  baseMean <- rowMeans(normCountsTable)
  
  ## baseMeanA
  baseMeanA <- rowMeans(normCountsTable[,groupA_Names])
  
  ## baseMeanB
  baseMeanB <- rowMeans(normCountsTable[,groupB_Names])
  
  ## foldChange
  foldChange <- baseMeanB/baseMeanA
  
  ## log2FoldChange
  log2FoldChange <- log2(foldChange)
  
  ## annotation
  # annotate BSU number 
  productDescr <- vector()
  for (BSU in geneIDs) {
    productDescr <- c(productDescr, annotTable[annotTable$BSU==BSU, ]$productVector)
  }
  
  DEGsTable <- data.frame(GeneID= geneIDs, productDescr=productDescr, baseMean= baseMean, baseMeanA=baseMeanA, 
                          baseMeanB=baseMeanB, foldChange=foldChange, log2FoldChange=log2FoldChange, 
                          pval= pvalsGLM, padj=padjGLM)
  return(DEGsTable)
}
# 
# ### Build an annotation (definition) vector for each BSU in the GeneID column
# definition <- vector()
# for (geneID in geneIDs) {
#   definition <- c(definition, keggGet(paste0("bsu:" ,geneID))[[1]]$DEFINITION)
# }
# 
# # the last enstry in geneIDs is the puratos gene that will not be found
# # Add one manually
# definition <- c(definition, "Puratos")

######################## $$$$$$$$$$$$$$$$$$$ ###################
outPath = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/DESeqOne/tables/Block4-Set2/"

############# 
# Contrast pair 1: 3 and 24 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==3]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==24]

DEGsTable_3hvs24h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_3hvs24h, paste0(outPath, "DESeqTable_3h.VS.24h.txt"), row.names = F, quote = F, sep="\t")

############# 
# Contrast pair 2: 3 and 48 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==3]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==48]

DEGsTable_3hvs48h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_3hvs48h, paste0(outPath, "DESeqTable_3h.VS.48h.txt"), row.names = F, quote = F, sep="\t")

############# 
# Contrast pair 3: 3 and 65 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==3]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==65]

DEGsTable_3hvs65h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_3hvs65h, paste0(outPath, "DESeqTable_3h.VS.65h.txt"), row.names = F, quote = F, sep="\t")

############# 
# Contrast pair 4: 24 and 48 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==24]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==48]

DEGsTable_24hvs48h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_24hvs48h, paste0(outPath, "DESeqTable_24h.VS.48h.txt"), row.names = F, quote = F, sep="\t")

############# 
# Contrast pair 5: 24 and 65 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==24]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==65]

DEGsTable_24hvs65h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_24hvs65h, paste0(outPath, "DESeqTable_24h.VS.65h.txt"), row.names = F, quote = F, sep="\t")

############# 
# Contrast pair 6: 48 and 65 hours
groupA_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==48]
groupB_Names <- blockFourDesignFile$Sample[blockFourDesignFile$Time==65]

DEGsTable_48hvs65h <- CalculateTable(cds2, groupA_Names, groupB_Names)
write.table(DEGsTable_48hvs65h, paste0(outPath, "DESeqTable_48h.VS.65h.txt"), row.names = F, quote = F, sep="\t")


################## $$$$$$$$$$$$$$$$$$$$$ ####################
##### $$$$$ #####
# Compile a list with all BSUs per GO term

BsubtilisGOA = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/6.B_subtilis_168.goa.txt"
GOATable <- read.delim(BsubtilisGOA, header = F, comment.char = "#")

# Build BSUs list 
# First Separate on pipe ("|"), grep only the elements containing BSUs and separate again on "/"
BSUs <- sapply(strsplit(as.character(GOATable$V11), "|", fixed=T), function(x) {
  index <- grep("BSU", x, fixed = T)
  #x[index]
  paste0("bsu:", unlist(strsplit(x[index], "/", fixed=T)))
})

allBSUsInGO <- unique(unlist(BSUs))


# extract all go terms present in bacillus annotation table and asign them to the BSUs list. 
# This way a GO => BSU mapping is build out of the annotation table
GOTerms <- GOATable$V5
names(BSUs) <- GOTerms

# aggregate all BSU that map to a GO term together so we build a GO => [BSUs] 
allGoTerms <- unique(GOTerms)
GOToBSU <- NULL
i = 1
for (term in allGoTerms){
  BSUVector <- unique(unlist(BSUs[ names(BSUs) == term ]))
  GOToBSU[[i]] <- BSUVector
  i = i + 1
}
names(GOToBSU) <- allGoTerms

# keep in pathways with more than 5 genes and less than 500
filterGOToBSU <- unlist(lapply(GOToBSU, function(x){
  if (length(x)>5 & length(x)<500) {
    T
  }else {
    F}
}))
GOToBSUFiltered <- GOToBSU[filterGOToBSU]

allBSUsInGOFiltered <- unique(unlist(GOToBSUFiltered))
length(allBSUsInGOFiltered)
length(allBSUsInGO)

######################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ##############################
# Pathway analysis. Use the code from Paul. 
##### Kegg pathways
# R-3.0.1
# Use hypergeometric test 

#biocLite("KEGGREST")
library("KEGGREST")
library("multtest")

listDatabases()
ls("package:KEGGREST") 


# bsu genes
bsu.genes <- keggList("bsu")
bsu.genes <- names(bsu.genes)
length(unique(bsu.genes)) # [1] 4421

# combine genes from Kegg and GO into one vector
BSUsKeggAndGO <- unique(c(bsu.genes, allBSUsInGO))

# bsu pathways
bsu.pw <- keggList("pathway", "bsu")
length(bsu.pw) # [1] 114

# Get genes per pathway
genes.pw <- list()
for(i in 1:length(bsu.pw)){
  temp <- keggGet(names(bsu.pw[i]))
  genes.temp <-  temp[[1]]$GENE
  if(length(genes.temp) > 0){
    genes.temp <- genes.temp[seq(1, length(genes.temp), by=2)]
    genes.pw[[names(bsu.pw[i])]] <- paste("bsu:", genes.temp, sep="")
    #genes.pw[[names(bsu.pw[i])]] <- genes.temp
  }
}
length(genes.pw) # [1] 102

#  

KeggAndGO_ToBSUs <- c(genes.pw, GOToBSUFiltered)

DEGs <- row.names(DEGsTable_3hvs24h[ DEGsTable_3hvs24h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 

# Over representation
hyper.test <- function(gene.set, selected, gene.list)
  # gene.set = pathway genes
  # selected = DEGs
  # gene.list = Background
{
  numWdrawn <- length(intersect(gene.set, selected))
  numDrawn <- length(selected)
  numW <- length(gene.set)
  numB <- length(gene.list) - numW
  
  #over represented
  pvalue <- phyper(numWdrawn - 1, numW, numB, numDrawn, lower.tail=FALSE) #P(X > x)
  
  #under represented
  #pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=TRUE) #P(X <= x)
  
  x <- max(0, numDrawn-numB):min(numW, numDrawn)
  y <- dhyper(x, numW, numB, numDrawn)
  names(y) <- x
  #return(list(observed=numWdrawn, estimated=y, pvalue=pvalue))
  return(pvalue)
} 


pw <- KeggAndGO_ToBSUs

####### 3h vs 24h
DEGs <- row.names(DEGsTable_3hvs24h[ DEGsTable_3hvs24h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <- data.frame(row.names = names(Overrep)) 
a <-cbind(a, adjPvals[ ,2])
colnames(a) <- "3hvs24h"

####### 3h vs 48h
DEGs <- row.names(DEGsTable_3hvs48h[ DEGsTable_3hvs48h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <-cbind(a, adjPvals[ ,2])
colnames(a)[2] <- "3hvs48h"

####### 3h vs 65h
DEGs <- row.names(DEGsTable_3hvs65h[ DEGsTable_3hvs65h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <-cbind(a, adjPvals[ ,2])
colnames(a)[3] <- "3hvs65h"

####### 24h vs 48h
DEGs <- row.names(DEGsTable_24hvs48h[ DEGsTable_24hvs48h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <-cbind(a, adjPvals[ ,2])
colnames(a)[4] <- "24hvs48h"

####### 24h vs 65h
DEGs <- row.names(DEGsTable_24hvs65h[ DEGsTable_24hvs65h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <-cbind(a, adjPvals[ ,2])
colnames(a)[5] <- "24hvs65h"

####### 48h vs 65h
DEGs <- row.names(DEGsTable_48hvs65h[ DEGsTable_48hvs65h$padj< 0.05, ])
DEG.names <- paste0("bsu:", DEGs) 
Overrep <- unlist(lapply(pw,hyper.test,selected=as.vector(DEG.names), gene.list=BSUsKeggAndGO ))
pvals <- as.vector(Overrep)
T.res <- mt.rawp2adjp(pvals, "BH")
adjPvals <- T.res$adjp[order(T.res$index),]

a <-cbind(a, adjPvals[ ,2])
colnames(a)[6] <- "48hvs65h"

#write.table(a, file=paste0(pathDesign, "pathwayAnalysis.txt"), sep="\t")


#################### $$$$$$$$$$$$$$ ####################
# annotate the GO terms 
library("GO.db")
library("gplots")

# a list of GO terms definitions (GO Term => Definition)
GODefinitions <- as.list(Definition(GOTERM))

# Combine the GO and Kegg terms Definition into one list
KeggAndGODefin <- c(bsu.pw, GODefinitions)

GOAndKeggDefinition <- vector()
GOAndKeggIDs <- row.names(a)
for (id in GOAndKeggIDs) {
  GOAndKeggDefinition <- c(GOAndKeggDefinition, KeggAndGODefin[id])
}

out <- data.frame(IDs = GOAndKeggIDs, Definitions= unlist(GOAndKeggDefinition), a)

write.table(out, file=paste0(pathDesign, "pathwayAnalysis.txt"), sep="\t", row.names = F)

#################### $$$$$$$$$$$$$$ ####################
# Heatmap for CCM pathway genes 

#genesGlycolysis <- genes.pw[["path:bsu00010"]]

# Fish out all "central carbihydrate metabolsim" genes
# First build a Module ID => Gene ID list
bsuModules <- keggList("module", "bsu")
length(bsuModules) # [1] 133

# Get genes per pathway
module2Genes <- list()
for (i in 1:length(bsuModules)) {
  temp <- keggGet(names(bsuModules[i]))
  genes.temp <-  temp[[1]]$GENE
  if (length(genes.temp) > 0){
    genes.temp <- genes.temp[seq(1, length(genes.temp), by=2)]
    module2Genes[[names(bsuModules[i])]] <- paste("bsu:", genes.temp, sep="")
    #genes.pw[[names(bsu.pw[i])]] <- genes.temp
  }
}
length(module2Genes) # [1] 102

# Compile all genes in CCM 
# List of Module IDs in CCM:
modulesCCM <- paste0("md:bsu_", c("M00001", "M00002", "M00307", "M00009", "M00010", "M00011", "M00004", "M00006", "M00007",
                "M00580", "M00005", "M00008", "M00308", "M00633", "M00309"))

genesInCCM <- vector()
for (id in modulesCCM){
  geneIDs <- module2Genes[[id]]
  genesInCCM <- c(genesInCCM, geneIDs)
}

genesInCCM <- unique(genesInCCM)

# Normalize count table
sizeFactors.mad <- function (counts, locfunc = median){
  loggeomeans <- rowMeans(log(counts))
  apply(counts, 2, function(cnts) exp(locfunc((log(cnts) -
                                                 loggeomeans)[is.finite(loggeomeans)])))
}
sf <- sizeFactors.mad(blockFourCountTable)

#divide countdata by sizefactors#

blockFourCountTable.scaled <- blockFourCountTable
for(i in 1:ncol(blockFourCountTable.scaled)){
  blockFourCountTable.scaled[,i] <- blockFourCountTable.scaled[,i]/sf[i]
} 
# Log2 the table
logBlockFourCountTable.scaled <- log2(blockFourCountTable.scaled)

## calculate z-scores from count table
zScoreCountTable = t(scale(t(blockFourCountTable.scaled))) 
head(zScoreCountTable)


row.names(blockFourCountTable) <- paste0("bsu:", row.names(blockFourCountTable))
glycolysisCountTable <- zScoreCountTable[paste0("bsu:", row.names(zScoreCountTable)) %in% genesInCCM, ]

# substitute BSUs with gene names in the glycolysis count table
BSUs <- paste0("bsu:", row.names(glycolysisCountTable))
geneNames <- vector()
for (BSU in BSUs) {
  geneNames <- c(geneNames, keggGet(BSU)[[1]]$NAME)
}
row.names(glycolysisCountTable) <- geneNames

timeFermentor <- paste0(blockFourDesignFile$Time, "h_", blockFourDesignFile$Culture)
colnames(glycolysisCountTable) <- timeFermentor

cols = c("#006400", "#0B6A00", "#157100", "#207700", "#2A7E00", "#358400", "#408B00",
          "#4A9100", "#559800", "#609E00", "#6AA500", "#75AB00", "#80B100", "#8AB800",
          "#95BE00", "#9FC500", "#AACB00", "#B5D200", "#BFD800", "#CADF00", "#D5E500",
          "#DFEC00", "#EAF200", "#F4F900", "#FFFF00", "#FFFF00", "#FFF400", "#FFEA00",
          "#FFDF00", "#FFD500", "#FFCA00", "#FFBF00", "#FFB500", "#FFAA00", "#FF9F00",
          "#FF9500", "#FF8A00", "#FF8000", "#FF7500", "#FF6A00", "#FF6000", "#FF5500",
          "#FF4A00", "#FF4000", "#FF3500", "#FF2A00", "#FF2000", "#FF1500", "#FF0B00",
          "#FF0000")

br = c(-3,seq(-2.5,2.5,length.out=49),3) 

# Draw the Heatmap
heatmap.2(glycolysisCountTable, dendrogram="row" , breaks=br,  key = T, col = cols, trace="none", Colv=NA)

### Build a table with all genes names and definitions in the CCM, as used for this analysis.
geneDefinitions <- vector()
for (geneID in genesInCCM){
  geneDefinitions <- c(geneDefinitions, keggGet(geneID)[[1]]$DEFINITION)
}
genesTableCCM <- data.frame(Name = geneNames, Definition = geneDefinitions)
write.table(genesTableCCM, paste0(outPath, "CCMGenesTable.txt"), row.names = F, sep = "\t", quote = F)

#adjPvals <- T.res$adjp[order(T.res$index),]
#rownames(adjPvals) <- names(Overrep)
#Overrep <- names(Overrep[which(Overrep < 0.1)])
#Overrep_adj <- adjPvals[as.vector(which(adjPvals[,2] < 0.1)),]


