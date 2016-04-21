# Author: Ioannis Moustakas, i.moustakas@uva.nl 
# Title: Combine count tables for puratos 

path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/mappingWithTophatNewData/MAPQTen/"

listOfFiles = grep("countMinQual*", dir(path), value=TRUE)

samples <- read.delim(paste0(path, listOfFiles[1]), header=F)

for (file in listOfFiles[2:length(listOfFiles)]){
  sample <- read.delim(paste0(path, file), header=F)
  samples <- cbind(samples, sample[,2])
}

colnames(samples) <- c("GeneName", listOfFiles)
write.table(samples, file=paste0(path, "allSamplesCountTable.txt"), sep="\t", row.names = F)

############################## $$$$$$$$$$$$$$$$$$$$$$$$ #################################
path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Scratch/alignOnGenome/"
readLength <- read.delim(paste0(path, "lengthReadsrrnB-23S.txt"), header=F)

readsLength <- readLength[,1]

hist(readsLength, breaks = 100)