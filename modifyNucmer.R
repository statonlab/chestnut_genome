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

## merge the combine maps to Allmaps
allmap <- read.csv("AllMaps_CM_May2018v6.csv", header = T)
combineQl_MA$CM <- gsub("lcl_","",combineQl_MA$CM)
names(combineQl_MA)[1] <- "Chestnut_contig"
new_allmap <- join(allmap[,-c(27,28)], combineQl_MA[,-8], by = "Chestnut_contig", match="first")
write.csv(new_allmap, "AllMaps_CM_Sep2019.csv")

QL_MA_only <- combineQl_MA[which(!combineQl_MA$Chestnut_contig %in% intersect(allmap$Chestnut_contig, combineQl_MA$Chestnut_contig)),]
write.csv(QL_MA_only, "contigs_not_in_allmap.csv")

## compare Ql and mahogany to the version 4.1 arrangement
CmVersion4 <- read.csv("Castanea_mollissima_scaffolds_v4.1_contigPositions.csv", header = T)
length(unique(CmVersion4$cc_lcl_contig))
## Ql uniquely mapped contigs vs v4.1
rows_dup <- which(Ql_contig_len$CM %in% unique(Ql_contig_len$CM[duplicated(Ql_contig_len$CM)]),)
uniq_Ql_Cm <- Ql_contig_len[-rows_dup,]
V4_shared_Ql <- intersect(uniq_Ql_Cm$CM, CmVersion4$cc_lcl_contig)
sum(uniq_Ql_Cm$length[which(uniq_Ql_Cm$CM %in% V4_shared_Ql)])
# Mahogany uniquely mapped contigs vs v4.1
rows_dup <- which(Mah_contig_len$CM %in% unique(Mah_contig_len$CM[duplicated(Mah_contig_len$CM)]),)
uniq_Ma_Cm <- Mah_contig_len[-rows_dup,]
V4_shared_Mah <- intersect(uniq_Ma_Cm$CM, CmVersion4$cc_lcl_contig)
sum(uniq_Ma_Cm$length[which(uniq_Ma_Cm$CM %in% V4_shared_Mah)])

## add version 4.1 map to allmaps
CmVersion4$cc_lcl_contig <- gsub("lcl_","",CmVersion4$cc_lcl_contig)
names(CmVersion4)[1] <- "Chestnut_contig"
new_allmapV2 <- join(new_allmap, CmVersion4, by = "Chestnut_contig", match = "first")
write.csv(new_allmapV2, "AllMaps_CM_Sep2019v2.csv")

## split HB2, JB1, NK4 parants LG
# dedup allmaps
dedupMap <- unique(allmap[,1:14])
# split HB2
HB2_A1map <- unique(allmap[grep("1",allmap$HB2_LG),][,c(1,15,16,17)])
HB2_A2map <- unique(allmap[grep("2", allmap$HB2_LG),][,c(1,15,16,17)])
new_allmapV3 <- join(dedupMap,HB2_A1map, by = "Chestnut_contig", match = "all")
names(new_allmapV3)[15:17] <- c("HB2_P1_Marker","HB2_P1_LG", "HB2_P1_cM")
new_allmapV3 <- join(new_allmapV3,HB2_A2map, by = "Chestnut_contig", match = "all")
names(new_allmapV3)[18:20] <- c("HB2_P2_Marker","HB2_P2_LG", "HB2_P2_cM")

# split JB1
JB1_A1map <- unique(allmap[grep("[A-Z]_", allmap$JB1_Marker),][,c(1,18,19,20)])
JB1_A2map <- unique(allmap[grep("[A-Z]2", allmap$JB1_LG),][,c(1,18,19,20)])
new_allmapV3 <- join(new_allmapV3, JB1_A1map, by = "Chestnut_contig", match = "all")
names(new_allmapV3)[21:23] <- c("JB1_P1_Marker","JB1_P1_LG","JB1_P1_cM")
new_allmapV3 <- join(new_allmapV3, JB1_A2map, by="Chestnut_contig", match= "all")
names(new_allmapV3)[24:26] <- c("JB1_P2_Marker","JB1_P2_LG","JB1_P2_cM")

# split NK4
NK4_A1map <- unique(allmap[grep("[A-Z]_", allmap$NK4_Marker),][,c(1,21,22,23)])
NK4_A2map <- unique(allmap[grep("2", allmap$NK4_LG),][,c(1,21,22,23)])
new_allmapV3 <- join(new_allmapV3, NK4_A1map, by="Chestnut_contig",match="all")
names(new_allmapV3)[27:29] <- c("NK4_P1_Marker","NK4_P1_LG","NK4_P1_cM")
new_allmapV3 <- join(new_allmapV3, NK4_A2map, by="Chestnut_contig",match = "all")
names(new_allmapV3)[30:32] <- c("NK4_P2_Marker","NK4_P2_LG","NK4_P2_cM")

# add the rest of mapps(oak_markers, v4.1_LG, Mahogany, Qlobata)
oak_map <- unique(allmap[,c(1,24,25,26)])
new_allmapV3 <- join(new_allmapV3, oak_map, by = "Chestnut_contig", match = "all")
new_allmapV3 <- join(new_allmapV3, CmVersion4, by = "Chestnut_contig", match="all")
new_allmapV3_wQM <- join(new_allmapV3, combineQl_MA[,-8], by = "Chestnut_contig", match="first")

write.csv(new_allmapV3_wQM, "AllMaps_CM_Sep2019v3.csv")
