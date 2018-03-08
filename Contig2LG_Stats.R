### Task
## To take the current set of contig to LG mappings
## And print out
## the number of contigs that map to a single LG across all maps
## their total length
## the number of contigs (and a list) that map to more than one LG 
## their total length
library(reshape2)

setwd("/Users/mestato/Desktop/chestnut_genome")

# read in data frame of all contig mappings
contig2lg <- read.csv("AllMaps_20180306.csv", stringsAsFactors = FALSE)

# remove columns we don't need for this
contig2lg <- contig2lg[,c("Chestnut_contig","Fan_LG","Kub_LG","JB1_LG","NK4_LG","HB2_LG")]

# remove rows with NAs for all LG columns
contig2lg <- contig2lg[rowSums(is.na(contig2lg)) < 5, ]

# remove numbers and whitespace from linkage groups
contig2lg$JB1_LG <- gsub("\\d+","", contig2lg$JB1_LG)
contig2lg$NK4_LG <- gsub("\\d+","", contig2lg$NK4_LG)
contig2lg$HB2_LG <- gsub("\\d+","", contig2lg$HB2_LG)

contig2lg$Fan_LG <- gsub("\\s+","", contig2lg$Fan_LG)
contig2lg$Kub_LG <- gsub("\\s+","", contig2lg$Kub_LG)
contig2lg$JB1_LG <- gsub("\\s+","", contig2lg$JB1_LG)
contig2lg$NK4_LG <- gsub("\\s+","", contig2lg$NK4_LG)
contig2lg$HB2_LG <- gsub("\\s+","", contig2lg$HB2_LG)

# remove duplicate rows
contig2lg <- unique(contig2lg)

# melt to merge LG columns into one giant column
contig2lg <- melt(contig2lg, id=(c("Chestnut_contig")))

# remove variable column - we don't care which map was used
contig2lg <- contig2lg[,c("Chestnut_contig","value")]

# remove rows with NA
contig2lg <- na.omit(contig2lg)

# remove duplicate rows
contig2lg <- unique(contig2lg)

# sort by contig name and renumber rows
contig2lg <- contig2lg[order(contig2lg$Chestnut_contig),]
row.names(contig2lg) <- 1:nrow(contig2lg)

