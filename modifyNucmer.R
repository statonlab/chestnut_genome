setwd("~/Desktop/Jiali/UTK/chestnut/blast_oak/")
# remove duplicate pairs
result <- read.table("Cm2Ql_CDS_filtered_annot_2708.txt", header = F)
names(result) <- c("CmID","Cm_contig","Cm_start","Cm_end","QlID","Ql_chr","Ql_start","Ql_end","BLASTidentity")

filter_result <- unique(result[,-9])
write.csv(filter_result, "filtered_Cm_Ql_geneBLAST.csv")
# Write the unique result into a csv file
length(unique(filter_result$CmID))

## nucmer chestnut v3.2 against mahogany and valley oak
library(stringr)
library(plyr)
nucmerQl_CM <- read.table("./Ql_Cm.delta.coords.filter",header = F,sep = "|", skip = 5)
nucmerMah_CM <- read.table("./Mahogany_Cm.coords.filter",header = F,sep = "|", skip = 5)
colnames(nucmerQl_CM) <- c("ref","query","length","identity","contig length","coverage","Tags")
colnames(nucmerMah_CM) <- c("ref","query","length","identity","contig length","coverage","Tags")

nucmerQl_CM$ref_end <- as.numeric(word(nucmerQl_CM$ref, -3)) 
nucmerQl_CM$AL <- as.numeric(word(nucmerQl_CM$length, -3))
nucmerQl_CM$ref_start <- nucmerQl_CM$ref_end - nucmerQl_CM$AL

## Filter contigs mapped to valley oak chromosomes
pairs <- unique(nucmerQl_CM$Tags)
filteredQl_Cm <- data.frame(matrix(ncol = 1, nrow = 0))
for (pair in pairs) {
rows_pair <- nucmerQl_CM[which(nucmerQl_CM$Tags == pair),]
rows_pair$query_len <- as.numeric(word(rows_pair$`contig length`,-3))
rows_pair$AL <- as.numeric(word(rows_pair$length,-3))
CALP <- sum(rows_pair$AL/rows_pair$query_len) *100
if (CALP > 5) {
filteredQl_Cm <- rbind(filteredQl_Cm, as.data.frame(pair))  
}}
colnames(filteredQl_Cm) <- c("Tags")
## Add start and end positions
filteredQl_Cm_coor <- join(filteredQl_Cm, nucmerQl_CM[order(nucmerQl_CM$ref_start), c("Tags", "ref_start")], by="Tags", match="first")
filteredQl_Cm_coor <- join(filteredQl_Cm_coor, nucmerQl_CM[order(nucmerQl_CM$ref_end, decreasing = T), c("Tags", "ref_end")], by="Tags", match="first")

# split the tags into chestnut contigs and oak
filteredQl_Cm_coor$Ql <- t(as.data.frame(str_split(filteredQl_Cm_coor$Tags,'\t')))[,1]
filteredQl_Cm_coor$CM <- t(as.data.frame(str_split(filteredQl_Cm_coor$Tags,'\t')))[,2]
filteredQl_Cm_coor$Ql_chr <- word(filteredQl_Cm_coor$Ql, -1)
filteredQl_Cm_coor <- filteredQl_Cm_coor[,c("CM","Ql_chr","ref_start","ref_end")]
length(unique(filteredQl_Cm_coor$CM))

## Filter contigs mapped to mahogny chromosomes
pairs <- unique(nucmerMah_CM$Tags)
filteredMah_Cm <- data.frame(matrix(ncol = 1, nrow = 0))
for (pair in pairs) {
  rows_pair <- nucmerMah_CM[which(nucmerMah_CM$Tags == pair),]
  rows_pair$query_len <- as.numeric(word(rows_pair$`contig length`,-3))
  rows_pair$AL <- as.numeric(word(rows_pair$length,-3))
  CALP <- sum(rows_pair$AL/rows_pair$query_len) *100
  if (CALP > 5) {
    filteredMah_Cm <- rbind(filteredMah_Cm, as.data.frame(pair))  
  }}
colnames(filteredMah_Cm) <- c("Tags")
## Add start and end positions
nucmerMah_CM$Ma_end <- as.numeric(word(nucmerMah_CM$ref, -3)) 
nucmerMah_CM$AL <- as.numeric(word(nucmerMah_CM$length, -3))
nucmerMah_CM$Ma_start <- nucmerMah_CM$Ma_end - nucmerMah_CM$AL

filteredMah_Cm_coor <- join(filteredMah_Cm, nucmerMah_CM[order(nucmerMah_CM$Ma_start), c("Tags", "Ma_start")], by="Tags", match="first")
filteredMah_Cm_coor <- join(filteredMah_Cm_coor, nucmerMah_CM[order(nucmerMah_CM$Ma_end, decreasing = T), c("Tags", "Ma_end")], by="Tags", match="first")

# split the tags into chestnut contigs and oak
filteredMah_Cm_coor$Ma <- t(as.data.frame(str_split(filteredMah_Cm_coor$Tags,'\t')))[,1]
filteredMah_Cm_coor$CM <- t(as.data.frame(str_split(filteredMah_Cm_coor$Tags,'\t')))[,2]
filteredMah_Cm_coor$Ma_chr <- word(filteredMah_Cm_coor$Ma, -1)
filteredMah_Cm_coor <- filteredMah_Cm_coor[,c("CM","Ma_chr","Ma_start","Ma_end")]
length(unique(filteredMah_Cm_coor$CM))

combineQl_MA <- merge(filteredMah_Cm_coor, filteredQl_Cm_coor, by="CM", all=T)
write.csv(combineQl_MA, "Cm_Ma_Ql_nucmer.csv")

# get seq len of the contiges
contig_len <- read.table("Castanea_mollissima_scaffolds_v3.2.lens.tsv", header = F)
seq_len <- sum(contig_len$V2)
seq_len
names(contig_len) <- c("CM", "length")
contig_len$CM <- gsub("\\|","_", contig_len$CM)
combineQl_MA <- join(combineQl_MA, contig_len, by = "CM", type="left")
Ql_contig_len <- join(filteredQl_Cm_coor, contig_len, by = "CM", type = "left")
Mah_contig_len <- join(filteredMah_Cm_coor, contig_len, by = "CM", type = "left")

# Calculate the number of contigs and length of mapped contigs   
stats_contigs <- function(df){
rows_dup <- which(df$CM %in% unique(df$CM[duplicated(df$CM)]),)
length(rows_dup)
length(df$CM) - length(rows_dup)
multimapped <- sum(df$length[rows_dup])
print(length(unique(df$CM[rows_dup])))
print(sum(unique(df$length[rows_dup])))
unqiue_len <- sum(df$length) - multimapped
print(length(df$CM) - length(df$CM[rows_dup]))
print(unqiue_len)
}
stats_contigs(Ql_contig_len)
stats_contigs(Mah_contig_len)
