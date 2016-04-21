# Title: annotate Puratos RNAseq count table. Should work for any other project
# Prerequizites: Counting is done on the basis of BSU IDs
# Author: Ioannis Moustakas, i.moustakas@uva.nl

projectPath='/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Scratch/alignOnGenomeNewDtata/'
gffFile='BSubltilis168_NCBIPlusPuratos.gff'
options(stringsAsFactors = FALSE)

# load the gff file 
gffTable <- read.delim(paste0(projectPath, gffFile), comment.char="#", stringsAsFactors=F)
colnames(gffTable) <- c("Contig", "RefSeq", "region", "Start", "End", "X0", "Strand", "X0.1",   "Attributes")

# extract the gene entries and put them in a table
gffTableGene <- gffTable[ gffTable$region=="gene", ]

# extract gene start and end coordinates, and the gene strand
GeneStart <- gffTableGene$Start
GeneEnd <- gffTableGene$End
GeneStrand <- gffTableGene$Strand

# Build an attribute table
GeneAttributes <- gffTableGene$Attributes
splitAttribures <- strsplit(GeneAttributes,  split=';', fixed=T)
#firstElement <- strsplit(splitAttribures[[1]],  split='=', fixed=T)
#header <- sapply(firstElement, function(x) x[1])
#geneAttributesTable <- data.frame(matrix(header, ncol=length(header)), stringsAsFactors=F)
#geneAttributesTable <- geneAttributesTable[-1, ]

# save the atributes into a list of data.frames 
geneAttributesList <- list()


for (i in 1:length(splitAttribures)){
  row <- strsplit(splitAttribures[[i]],  split='=', fixed=T)
  attributeValues <- sapply(row, function(x) toString(x[2]))
  attributeNames <- sapply(row, function(x) toString(x[1]))
  element <- as.data.frame(matrix(attributeValues, ncol=length(attributeValues)))
  colnames(element) <- attributeNames
  geneAttributesList[[i]] <- element
 # geneAttributesTable <- rbind(geneAttributesTable, attributeValues)
}

# go through the list of data.frames and build the annotation table
geneAnnotationTable <- do.call(rbind, lapply(geneAttributesList, function(x) {
  data <- c(x$Name, x$locus_tag, x$ID)
  matrix(data, ncol=3)
}))

geneAnnotationDF <- as.data.frame(geneAnnotationTable)
anntationTable <- cbind(geneAnnotationDF, GeneStart, GeneEnd, GeneStrand)
colnames(anntationTable) <- c("GeneID", "BSU", "Parent", "GeneStart", "GeneEnd", "GeneStrand")

############### $$$$$$$$$$$$$$ ###############
# Now extract the Product description form the CDS fields
gffTableCDS <- gffTable[ gffTable$region=="CDS", ]
CDSAttributes <- gffTableCDS$Attributes
splitAttribures <- strsplit(CDSAttributes,  split=';', fixed=T)
firstElement <- strsplit(splitAttribures[[1]],  split='=', fixed=T)
header <- sapply(firstElement, function(x) x[1])
CDSattributesTable <- data.frame(matrix(header, ncol=length(header)), stringsAsFactors=F)
CDSattributesTable <- CDSattributesTable[-1, ]

for (i in 1:length(splitAttribures)){
  row <- strsplit(splitAttribures[[i]],  split='=', fixed=T)
  attributeValues <- sapply(row, function(x) toString(x[2]))
  CDSattributesTable <- rbind(CDSattributesTable, attributeValues)
}
colnames(CDSattributesTable) <- header

productAnnot <- data.frame(GeneID=CDSattributesTable$Parent, Product=CDSattributesTable$product)

BSUs <- anntationTable$Parent
productVector <- c()
for (id in BSUs){
  if (id %in%  productAnnot$GeneID ){
    product <- productAnnot[ productAnnot$GeneID == id, ]$Product
    productVector <- c(productVector, product[1])
  } else productVector <- c(productVector, NA)
}

annotationAppendix <- cbind(anntationTable, productVector)
  
write.table(annotationAppendix, file = paste0(projectPath, 'AnnotationAppendix.txt'), row.names = F, sep='\t' )
  
  
  
