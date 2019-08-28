setwd("~/Desktop/Jiali/UTK/chestnut/blast_oak/chestnut_genome/")
library(plyr)
filtered_Ql_blast <- read.table("Blast_Mummer/CmMarker2Ql_blastn_270819_filtered.txt", header = F)
names(filtered_Ql_blast) <- c('marker_name', 'sseqid', 'qlen','pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sstrand')

uniq_blast <- filtered_Ql_blast[!duplicated(filtered_Ql_blast[,1:2]),]
Kub <- read.csv("./maps_table/11295_2012_579_MOESM5_ESM.csv", header = T)
Kub$marker_name <- gsub(" ","",Kub$marker_name)
Ql_Kub <- join(Kub[,1:3], uniq_blast[,c(1,2,10,11)], by ="marker_name", type="left", match = "first")
write.csv(Ql_Kub, "Cm_marker_to_Ql.csv")

