# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: combine ERCCs read counts with gene read counts

# Path ewhere the ERCCs are stored
ERCCsPath="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophat/ERCCsAlignment/countERCCS"

# set the working path and make a list of all text files in it
listOfFiles <- grep(".txt$", dir(ERCCsPath), value=TRUE)
# read the first file in the list and put in the 
ERCCsTableFirst <- read.delim(listOfFiles[1], header=F)
ERCCsAllSamples <- data.frame(Names=ERCCsTableFirst[,1])

# Go through all the files in the folder, and add them to the Data Frame
for (file in listOfFiles) {
  ERCCsTable <- read.delim(file, header=F)
  # remove the last line, this is the unmapped reads
  ERCCsTable <- ERCCsTable[-nrow(ERCCsTable), ]
  ID <- strsplit(x = file, split = "_", fixed = T)[[1]][1]
  colnames(ERCCsTable) <- c("Names", "SeqLength", "ReadsMapped", "ReadsNotMapped")
  ERCCsNamesReads <- data.frame( ERCCsTable$Names, ERCCsTable$ReadsMapped)
  colnames(ERCCsNamesReads) <- c("Names", ID)
  ERCCsAllSamples <- merge(ERCCsAllSamples, ERCCsNamesReads, by="Names")
}

# Load the gene counts table
geneCounts="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophat/extraSeqMAPQTen/allSamplesCountTable.txt"
geneCountTable <- read.delim(geneCounts)

# repalce the header 
colnames(geneCountTable) <- colnames(ERCCsAllSamples)
ERCCsGeneCountTable <- rbind(ERCCsAllSamples, geneCountTable)
write.table(ERCCsGeneCountTable, 
            file="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophat/extraSeqMAPQTen/combinedERCCsGeneCountTable.txt", 
            sep="\t",
            row.names=F)

