# Title: Normalize the RNAseq count table
# Author: Ioannis Moustakas, i.moustakas@uva.nl, usning M. Jonker core code

# !!!!!!!!!!!!!!! first set the correct working dir !!!!!!!!!!!!!!!!
# load the table
countTable <- read.delim("allSamplesCountTableNoS08.csv")

# remove the bottom 5 lines that are not gene counts but the report of htseq-count
nrows <- nrow(countTable) 
countTableGenesOnly <- countTable[-c((nrows-4):nrows), ] 

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

# Write the normalized table to a file
write.table(CountTable.scaled, file='CountTableNormalizedNoS08.csv', sep='\t')
