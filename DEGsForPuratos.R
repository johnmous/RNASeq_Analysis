# Author: Ioannis Moustakas, i.moustakas@uva.nl (Based on Martijs Jonker DEGs Script, m.j.jonker@uva.nl )
# Title: Perform a DEG analysis
# Save the comparison pairs into a table and feed it as input. This is handy when a large number of comparisons is required.
# There is an annotation function that is serving to anontate the gene names with a short description of the product function  

setwd("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/")

options(stringsAsFactors = FALSE)

path <- c("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/DEG_Analysis/")

# load the design file 
design <- read.delim(paste0(path,"designDEGAnalysis.csv"),stringsAsFactors=FALSE)
# load the gene count table
df <- read.delim("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Results/mapping/allSamplesCountTable.txt",stringsAsFactors=FALSE)
# remove the bottom 5 lines that are not gene counts but the report of htseq-count
nrows <- nrow(df) 
ncolns <- ncol(df)
df <- df[-c((nrows-4):nrows), ] 
geneIDs <- as.vector(df[,1],mode="character")
CountTable <- data.matrix(df[,2:ncolns])
tmp <- colnames(CountTable)
#tmp <- as.vector(unlist(lapply(strsplit(tmp,"_"),function(x){x[1]})),mode="numeric")
design <- design[match(tmp,design$Sample),]

# Remove genes whose read count is ) across all samples
check <- rowSums(CountTable)
length(which(check==0))
out <- which(check==0)
CountTable <- CountTable[-out,]
geneIDs <- geneIDs[-out]



# png("../../Images/Ioannis/Boxplot.png",width=635,height=460)
#   #cols <- ifelse(design$sampleinfo=="Control","blue","red")
#   cols <- design$Group
#   boxplot(log2(CountTable+0.5),pch=".",names=design$time.point,las=2,col=cols)
# graphics.off()
# 
# png("../../Images/Ioannis/Densityplot.png",width=512,height=512)
#   x <- log2(CountTable+0.5)
#   cols <- design$Group
#   plot(density(x[,1]),type="n",ylim=c(0,0.25))
#   for(i in 1:ncol(x)){lines(density(x[,i]),col=cols[i])}
# graphics.off()
# 
# PCA <- prcomp(t(log2(CountTable+0.5)))
# colnames(PCA$x) <- paste(colnames(PCA$x), " (", round(summary(PCA)$importance[2,]*100,1), ")%",sep="")
# PCAdata <- cbind(design,PCA$x)
# write.table(PCAdata,"../../Images/PCA.raw.txt",sep="\t",row.names=FALSE,quote=FALSE)

################ $$$$$$$$$$$$$$$$ ##################
# Build an annnotation table
#  Gene Name => Note from GFF file
# Load the annotation file
annotation <- read.delim("/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/sly/DNA/genome/ITAG2.4_gene_models.gff3", header = F, comment.char = "#", stringsAsFactors = F)
# rename the columns of annotation table
names(annotation)[1] <- "seqname"
names(annotation)[2] <- "RefSeq"
names(annotation)[3] <- "feature"
names(annotation)[4] <- "start"
names(annotation)[5] <- "end"
names(annotation)[7] <- "strand"
names(annotation)[9] <- "attribute"

mRNAFields <- annotation[ annotation$feature=="mRNA", ]

# A function to extract the gene description from the attribute field of gff file
extractGeneAttribute <- function(attribute) {
  names <- c()
  notes <- c()
  for (attrLine in attribute){
    if (length(attrLine)==0){
      return(NA)
    } else { 
      nameField <- strsplit(attrLine, ";", fixed=T)[[1]]
      # key (attribbute tag) => value (attribute value), stored in a Data Frame
      df <- data.frame()
      for (attr in nameField){ 
        key <- strsplit(attr, "=", fixed=T)[[1]][1]
        value <- strsplit(attr, "=", fixed=T)[[1]][2] 
        line <- data.frame(value=value)
        row.names(line) <- key
        df <- rbind(df, line)
      }
      # 
      name <- as.character(df["Name", ])
      name <- substr(name, 1, nchar(name)-2)
      names <- c(names, name)
      notes <- c(notes, as.character(df["Note", ]))
    }
  }
  namesNotes <- data.frame(Names = names, Notes = notes)
  return(namesNotes)
}


name2NoteDf <- extractGeneAttribute(mRNAFields$attribute)

AnnotationVector <- function(res, name2NoteDf) {
  genesList <- strsplit(res$id, split="gene:", fixed=T)
  genes <- sapply(genesList, function(x) x[2])
  annotation <- sapply(genes, function(x) {
    note <- as.character(name2NoteDf[name2NoteDf$Names==x, ]$Notes)
    if (length(note)==0) note <- NA
    note
  })
  return(annotation)
}

######### $$$$$$$$$$ #########
# set working directories
tablePath = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/plotsAndTables/comparissonTables/"
plotPath = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/plotsAndTables/plots/"
dir.create(tablePath)

# Build a vector where the line of tomato plant and time point are combined
# this will be used for the selection of the experiment
condition <- factor(paste(design$line.of.tomato.plant, design$time.point, sep="_"))

# create data.frame to store some sample statistics
statistics <- data.frame(comparisson=c(1), pairContrast=c(1), minPadj=c(1), numGenesPadjBelThres=c(1))

# Load the table with the comparissons to make
allComparissons <- read.delim("./plotsAndTables/comparissonsMatrix.csv")

library(DESeq)
ls("package:DESeq")

rownames(CountTable) <- geneIDs
cds <- newCountDataSet(CountTable,condition)

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))

for (i in 1:nrow(allComparissons)) {
  x <- as.character(allComparissons[i, ])
  # idx <- which((condition==x[2]) | (condition==x[3]))
  fileName = paste(x[2], x[3], sep = ".VS.")
  
  # tmp <- cds[,idx]
  tmp <- cds
  tmp <- estimateDispersions(tmp)
  
  # create Dir
  dir.create( paste0(plotPath, x[1]), showWarnings = F)
  
  # plot  Dispertion estimates
  completePlotPath= paste0(plotPath, x[1],"/DispEstPlot_", fileName, ".png")
  png(completePlotPath, width=635,height=460)
  plotDispEsts(tmp)
  graphics.off()
  
  tmp@phenoData@data$condition
  res <- nbinomTest(tmp, x[2], x[3])
  # head(res)
  
  completePlotPath= paste0(plotPath, x[1],"/MAplot_", fileName, ".png")
  png(completePlotPath, width=635,height=460)
  plotMA(res)
  graphics.off()
  
  minPadj <- min(res$padj[!is.na(res$padj)])
  geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))
  comparissonStats <- data.frame(comparisson = x[1], pairContrast=fileName, minPadj=minPadj, numGenesPadjBelThres=geneCount)
  statistics <- rbind(statistics, comparissonStats)
  
  ## plot Padj Histogram
  
  completePlotPath <- paste0(plotPath, x[1],"/histPadj_", fileName, ".png")
  png(completePlotPath, width=635,height=460)
  hist(res$pval,breaks=30)
  graphics.off()
  
  # annotate Table
  annotation <- unlist(AnnotationVector(res, name2NoteDf))
  resAnnotated <- res
  resAnnotated$Annotation <- annotation
  resAnnotated <- resAnnotated[order(resAnnotated$padj), ]
  
  # create Dir
  dir.create( paste0(tablePath, x[1]), showWarnings = F)
  completeTablePath = paste0(tablePath, x[1],"/DESeqTable_", fileName, ".txt")
  write.table(resAnnotated, completeTablePath, sep="\t", row.names = F)
}

completeTablePath = paste0(tablePath, "cumulativeStatistics.txt")
write.table(statistics, completeTablePath, sep="\t", row.names = F)
