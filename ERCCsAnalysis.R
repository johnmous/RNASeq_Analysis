# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: Plot the ERCCs for Takken Experiment
library( ggplot2)
library( reshape2)
library(ggvis)
interactive()

# load the ERCCs concetration table
concetrationTable <- read.delim("/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/ERCC/ERCC_Controls_Analysis.txt", header=T)
ERCCsConcTable <- concetrationTable[,c(2,4)]
colnames(ERCCsConcTable)[1] <- "Names"

# set the working path and make a list of all text files in it
listOfFiles <- grep(".txt$", dir(), value=TRUE)
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



# put a column for the concetration of each of the ERCCs
ERCCsAllSamplesConc <- merge(ERCCsAllSamples, ERCCsConcTable, by="Names")
# Set the name of the ERCCs as the row name
row.names(ERCCsAllSamplesConc) <- ERCCsAllSamplesConc$Names
ERCCsAllSamplesConc[,1] <- ERCCsAllSamplesConc[,ncol(ERCCsAllSamplesConc)]
colnames(ERCCsAllSamplesConc)[1] <- "Concentration"
ERCCsAllSamplesConc <- ERCCsAllSamplesConc[,-ncol(ERCCsAllSamplesConc)]


ERCCsAllSamplesConcCollapsed <- ddply(ERCCsAllSamplesConc, "Concentration", numcolwise(sum))
test <- melt(ERCCsAllSamplesConcCollapsed, id="Concentration")
names(test) <- c("Concentration", "Sample", "Count")
# log trans
test$Count <- log2(test$Count+1)

add_tooltip(vis, html, on = c("hover", "click"))

all_values <- function(x) {
  if(is.null(x)) return(NULL)
  paste0(names(x), ": ", format(x), collapse = "<br />")
}

  
vis <- test %>% ggvis(~Concentration, ~Count, size.hover := 200) %>% 
  scale_numeric("x", trans="log", expand=0) %>% layer_points(fill=~factor(Sample)) %>% 
  add_tooltip(all_values, "hover")


base %>% add_tooltip(all_values, "click")




p <- ggvis(ERCCsAllSamplesConcCollapsed, x = ~log2(Concentration), y = ~log2(Count))
layer_points(p)
p <- ggvis(ERCCsAllSamplesConcCollapsed, x = ~log2(Concentration), y = ~log2(S02), stroke:= "red")
layer_points(p)

plot(PlotSpikeCounts(data, log_y_axis=T))

  