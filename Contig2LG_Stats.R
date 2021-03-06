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
contig2lg <- read.csv("AllMaps_20180501.csv", stringsAsFactors = FALSE)

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

###---------------------------------------------
### Subset df into contigs with a single LG (not chimeras) vs those with more than one LG (chimeras)
###---------------------------------------------

# get a logical vector with all duplicate contigs marked
contig2lg_duplicated <- duplicated(contig2lg$Chestnut_contig) | duplicated(contig2lg$Chestnut_contig, fromLast = TRUE)

# filter dataframe for those that are duplicated
contig2lg_chimera <- contig2lg[contig2lg_duplicated,]
contig_chimera <- contig2lg_chimera[,c("Chestnut_contig")]
contig_chimera <- unique(contig_chimera)


# and filter dataframe for those that are not duplicated
contig2lg_not_chimera <- contig2lg[!contig2lg_duplicated,]
contig_not_chimera <- contig2lg_not_chimera[,c("Chestnut_contig")]
contig_not_chimera <- unique(contig_not_chimera)

###---------------------------------------------
### Get contig length stats for each list (chimera and non chimera contigs)
###---------------------------------------------

# Read in contig lengths and merge with dfs
contig2len <- read.table("./maps_table/Castanea_mollissima_scaffolds_v3.2.lens.tsv", stringsAsFactors = FALSE,
                         col.names = c("Chestnut_contig", "Contig_length"))

# clean up contig names by removing lcl
contig2len$Chestnut_contig <- gsub("lcl\\|","", contig2len$Chestnut_contig)

# filter the len df by list of contigs
contig_chimera_lengths <- subset(contig2len, contig2len$Chestnut_contig %in% contig_chimera)
length(contig_chimera)
sum(contig_chimera_lengths$Contig_length)

contig_not_chimera_lengths <- subset(contig2len, contig2len$Chestnut_contig %in% contig_not_chimera)
length(contig_not_chimera)
sum(contig_not_chimera_lengths$Contig_length)

# print a summary of results
print( paste( 'Contigs that are mapped to a single LG: ',
               length(contig_not_chimera),
               '(',
              sum(contig_not_chimera_lengths$Contig_length),
               'bp)'
            )
     )

print( paste( 'Contigs that are chimeric (i.e. mapped to more than one LG): ',
              length(contig_chimera),
              '(',
              sum(contig_chimera_lengths$Contig_length),
              'bp)'
            )
      )

# create a file of chimeras - one file with the uniq list, one with the LGs
write.csv(contig_chimera_lengths, "contig_chimera_lengths_postChimeraRemoval0501.csv")
write.csv(contig2lg_chimera, "contig_chimera_LGs_0501_postChimeraRemoval0501.csv")

