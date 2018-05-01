setwd("/Users/mestato/Desktop/chestnut_genome")

library(plyr)
library(stringr)
##################### read maps in #####################

#----------------Fan's map----------------------------------
FanMap <- read.csv("./maps_table/Chinese_only_integrated_edit3_LG.txt", sep = '\t')
# convert to "CCall_contigXX_v2"
FanMap$Markers <- gsub("_\\d+","\\_v2", FanMap$Markers)
FanMap$Markers <- gsub("v2c","\\_contig", FanMap$Markers)
head(FanMap)

#-------------------------Kubisiak's map---------------------------
Kub <- read.csv("./maps_table/11295_2012_579_MOESM5_ESM.csv", header = T)
head(Kub)
colnames(Kub)[1] <- c("Markers")
Kub <- Kub[,1:3]
head(Kub)
Kub$Markers <- gsub("\\s+", "", Kub$Markers)

#--------------------HB2 map------------------------------
HB2 <- read.csv("./maps_table/HB2_map.csv", header = T, sep = '\t')
colnames(HB2)[3:4] <- c("name", "cM")
names(HB2)
HB2$HB2_Marker <- paste0(HB2$ID,"_",HB2$name)
names(HB2)

#---------------------JB1 map--------------------------------
JB1 <- read.csv("./maps_table/JB1map.csv", header=F)
JB1$JB1_Marker <- paste0(JB1$V2,"_",JB1$V4)

#--------------------NK4 map---------------------------------
NK4 <- read.csv("./maps_table/NK4map.csv", header=F)
NK4$NK4_Marker <- paste0(NK4$V2,"_",NK4$V4)



#--------------------------Oak map---------------------------
Oak <- read.csv("./maps_table/Composite_file_map_Vf_260413.csv", header = T)
Oak <- Oak[,1:3]
head(Oak)


#--------------------peach blast filtering-------------------------
blast_peach <- read.table("./Blast_Mummer_postChimeraRemoval/Pp2Cm_blastn.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_peach <- aggregate(V11~ V1+V2, data= blast_peach, FUN=min)
filterBLAST_peach <- filterBLAST_peach[order(filterBLAST_peach$V1),]
filterBLAST_peach$V1_V11 <- paste0(filterBLAST_peach$V1,".",filterBLAST_peach$V11)
rows_dup <- which(filterBLAST_peach$V1_V11 %in% unique(filterBLAST_peach$V1_V11[duplicated(filterBLAST_peach$V1_V11)]),)
dupdf <- filterBLAST_peach[rows_dup,]
filterBLAST_peach$V1 <- as.factor(filterBLAST_peach$V1)
noMultimatch_peach <- filterBLAST_peach[!filterBLAST_peach$V1 %in% as.character(dupdf$V1),]
noMultimatch_peach <- noMultimatch_peach[order(noMultimatch_peach$V11),]
peach_uniq <- as.data.frame(unique(noMultimatch_peach$V1))
colnames(peach_uniq) <- c("V1")
peach_uniq <- join(peach_uniq, noMultimatch_peach, by = "V1", match = "first")

## Add marker location on the contigs
peach_uniq$ID <- paste0(peach_uniq$V1,"_",peach_uniq$V2)
blast_peach$ID <- paste0(blast_peach$V1,"_",blast_peach$V2)
blast_peach <- blast_peach[order(blast_peach$V9),]
peach_uniq <- join(peach_uniq, blast_peach[,c("V9","ID")], by = "ID", match = "first")
blast_peach <- blast_peach[order(blast_peach$V10, decreasing = T),]
peach_uniq <- join(peach_uniq, blast_peach[,c("V10","ID")], by = "ID", match = "first")
names(peach_uniq)
peach_uniq <- join(peach_uniq, blast_peach[,c("V13", "ID")], by = "ID", match = "first")
colnames(peach_uniq)[6:7] <- c("start", "stop")
#write.csv(peach_uniq,"BLAST_Peach.csv")

##------------------------- Peach to chestnut mummer filtering--------------------------------
filterPp2_CM <- read.table("./Blast_Mummer_postChimeraRemoval/Pp2_Cm.coords.dir.filter", header = F, sep = "|", skip = 5)

colnames(filterPp2_CM) <- c("ref","query","length","identity","contig length","coverage", "Tags")

filterPp2_CM$Tmp <- t(as.data.frame(str_split(filterPp2_CM$Tags,'\t')))[,1]
filterPp2_CM$Tag_dir <- t(as.data.frame(str_split(filterPp2_CM$Tags,'[[:space:]]+')))[,3]
filterPp2_CM$Tag_peach <- t(as.data.frame(str_split(filterPp2_CM$Tags,'[[:space:]]+')))[,4]

filterPp2_CM$Tag_Cm <- t(as.data.frame(str_split(filterPp2_CM$Tags,'\t')))[,2]
peach2chestnut <- filterPp2_CM[,c("ref","Tags","Tag_peach","Tag_Cm", "Tag_dir")]
colnames(peach2chestnut) <- c("Location", "Tags", "Peach", "V2", "Dir")

## delete contigs mapped to multiple peach chromosomes
rows_dup_contigs <- which(filterPp2_CM$Tags %in% unique(filterPp2_CM$Tags[duplicated(filterPp2_CM$Tags)]),)
dupdf_contigs <- filterPp2_CM[rows_dup_contigs,]
uniq_pair <- filterPp2_CM[-c(rows_dup_contigs),]
deduped_contigs <- dupdf_contigs[-c(which(duplicated(dupdf_contigs$Tags),)),]
deduped_pair <- rbind(uniq_pair, deduped_contigs)

rows_multi_loc <- which(deduped_pair$Tag_Cm %in% unique(deduped_pair$Tag_Cm[duplicated(deduped_pair$Tag_Cm)]),)
uniq_contigs <- deduped_pair[-c(rows_multi_loc),]



#----------------Fan's map filtering--------------------------------
blast_fan <- read.table("./Blast_Mummer_postChimeraRemoval/markers_blastn.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_fan <- aggregate(V11~ V1+V2, data= blast_fan, FUN=min)
filterBLAST_fan <- filterBLAST_fan[order(filterBLAST_fan$V1),]
filterBLAST_fan$V1_V11 <- paste0(filterBLAST_fan$V1,".",filterBLAST_fan$V11)
rows_dup <- which(filterBLAST_fan$V1_V11 %in% unique(filterBLAST_fan$V1_V11[duplicated(filterBLAST_fan$V1_V11)]),)
dupdf <- filterBLAST_fan[rows_dup,]
filterBLAST_fan$V1 <- as.factor(filterBLAST_fan$V1)
noMultimatch_fan <- filterBLAST_fan[!filterBLAST_fan$V1 %in% as.character(dupdf$V1),]
noMultimatch_fan$V11_V2 <- paste0(noMultimatch_fan$V11,"|", noMultimatch_fan$V2)
noMultimatch_fan <- noMultimatch_fan[order(noMultimatch_fan$V11),]
Fan_uniq <- as.data.frame(unique(noMultimatch_fan$V1))
colnames(Fan_uniq) <- c("V1")
Fan_uniq <- join(Fan_uniq, noMultimatch_fan, by = "V1", match = "first")

## Add marker location on the contigs
Fan_uniq$ID <- paste0(Fan_uniq$V1,"_",Fan_uniq$V2)
blast_fan$ID <- paste0(blast_fan$V1,"_",blast_fan$V2)
blast_fan <- blast_fan[order(blast_fan$V9),]
Fan_uniq <- join(Fan_uniq, blast_fan[,c("V9","ID")], by = "ID", match = "first")
blast_fan <- blast_fan[order(blast_fan$V10, decreasing = T),]
Fan_uniq <- join(Fan_uniq, blast_fan[,c("V10","ID")], by = "ID", match = "first")
names(Fan_uniq)
colnames(Fan_uniq)[7:8] <- c("start", "stop")

#----------------HB2 map filtering------------------------------
blast_HB2 <- read.table("./Blast_Mummer_postChimeraRemoval/HB2_blastn.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_HB2 <- aggregate(V11~ V1+V2, data= blast_HB2, FUN=min)
filterBLAST_HB2 <- filterBLAST_HB2[order(filterBLAST_HB2$V1),]
filterBLAST_HB2$V1_V11 <- paste0(filterBLAST_HB2$V1,".",filterBLAST_HB2$V11)
rows_dup_H <- which(filterBLAST_HB2$V1_V11 %in% unique(filterBLAST_HB2$V1_V11[duplicated(filterBLAST_HB2$V1_V11)]),)
dupdf_H <- filterBLAST_HB2[rows_dup_H,]
filterBLAST_HB2$V1 <- as.factor(filterBLAST_HB2$V1)
noMultimatch_HB2 <- filterBLAST_HB2[!filterBLAST_HB2$V1 %in% as.character(dupdf_H$V1),]
noMultimatch_HB2$V11_V2 <- paste0(noMultimatch_HB2$V11,"|", noMultimatch_HB2$V2)
noMultimatch_HB2 <- noMultimatch_HB2[order(noMultimatch_HB2$V11),]
HB2_uniq <- as.data.frame(unique(noMultimatch_HB2$V1))
colnames(HB2_uniq) <- c("V1")
HB2_uniq <- join(HB2_uniq, noMultimatch_HB2, by = "V1", match = "first")

## Add marker location on the contigs
HB2_uniq$ID <- paste0(HB2_uniq$V1,"_",HB2_uniq$V2)
blast_HB2$ID <- paste0(blast_HB2$V1,"_",blast_HB2$V2)
blast_HB2 <- blast_HB2[order(blast_HB2$V9),]
HB2_uniq <- join(HB2_uniq, blast_HB2[,c("V9","ID")], by = "ID", match = "first")
blast_HB2 <- blast_HB2[order(blast_HB2$V10, decreasing = T),]
HB2_uniq <- join(HB2_uniq, blast_HB2[,c("V10","ID")], by = "ID", match = "first")
names(HB2_uniq)
colnames(HB2_uniq)[7:8] <- c("start", "stop")
colnames(HB2_uniq)[2]<- "HB2_Marker"

# Add LG and cM to the matched markers
HB2_new <- merge(HB2_uniq, HB2[,c(2,4,9)],by="HB2_Marker")
HB2_new <- HB2_new[,c(1,3,9,10)]
colnames(HB2_new)[2:4] <- c("contig","HB2_Group","HB2_cM")
HB2_new$HB2_Marker <- paste0("HB2_",HB2_new$HB2_Marker)

#------------------------JB1 map blast filtering--------------------------------------
blast_JB1 <- read.table("./Blast_Mummer_postChimeraRemoval/JB1_blastn.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_JB1 <- aggregate(V11~ V1+V2, data= blast_JB1, FUN=min)
filterBLAST_JB1 <- filterBLAST_JB1[order(filterBLAST_JB1$V1),]
filterBLAST_JB1$V1_V11 <- paste0(filterBLAST_JB1$V1,".",filterBLAST_JB1$V11)
rows_dup_H <- which(filterBLAST_JB1$V1_V11 %in% unique(filterBLAST_JB1$V1_V11[duplicated(filterBLAST_JB1$V1_V11)]),)
dupdf_H <- filterBLAST_JB1[rows_dup_H,]
filterBLAST_JB1$V1 <- as.factor(filterBLAST_JB1$V1)
noMultimatch_JB1 <- filterBLAST_JB1[!filterBLAST_JB1$V1 %in% as.character(dupdf_H$V1),]
noMultimatch_JB1$V11_V2 <- paste0(noMultimatch_JB1$V11,"|", noMultimatch_JB1$V2)
noMultimatch_JB1 <- noMultimatch_JB1[order(noMultimatch_JB1$V11),]
JB1_uniq <- as.data.frame(unique(noMultimatch_JB1$V1))
colnames(JB1_uniq) <- c("V1")
JB1_uniq <- join(JB1_uniq, noMultimatch_JB1, by = "V1", match = "first")

## Add marker location on the contigs
JB1_uniq$ID <- paste0(JB1_uniq$V1,"_",JB1_uniq$V2)
blast_JB1$ID <- paste0(blast_JB1$V1,"_",blast_JB1$V2)
blast_JB1 <- blast_JB1[order(blast_JB1$V9),]
JB1_uniq <- join(JB1_uniq, blast_JB1[,c("V9","ID")], by = "ID", match = "first")
blast_JB1 <- blast_JB1[order(blast_JB1$V10, decreasing = T),]
JB1_uniq <- join(JB1_uniq, blast_JB1[,c("V10","ID")], by = "ID", match = "first")
names(JB1_uniq)
colnames(JB1_uniq)[7:8] <- c("start", "stop")
colnames(JB1_uniq)[2]<- "JB1_Marker"

# Add LG and cM to the matched markers
JB1_new <- merge(JB1_uniq, JB1[,c(2,5,10)],by="JB1_Marker")
JB1_new <- JB1_new[,c(1,3,9,10)]
colnames(JB1_new)[2:4] <- c("contig","JB1_Group","JB1_cM")

#----------------------NK4 map filtering ------------------------------------
blast_NK4 <- read.table("./Blast_Mummer_postChimeraRemoval/NK4_blastn.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_NK4 <- aggregate(V11~ V1+V2, data= blast_NK4, FUN=min)
filterBLAST_NK4 <- filterBLAST_NK4[order(filterBLAST_NK4$V1),]
filterBLAST_NK4$V1_V11 <- paste0(filterBLAST_NK4$V1,".",filterBLAST_NK4$V11)
rows_dup_H <- which(filterBLAST_NK4$V1_V11 %in% unique(filterBLAST_NK4$V1_V11[duplicated(filterBLAST_NK4$V1_V11)]),)
dupdf_H <- filterBLAST_NK4[rows_dup_H,]
filterBLAST_NK4$V1 <- as.factor(filterBLAST_NK4$V1)
noMultimatch_NK4 <- filterBLAST_NK4[!filterBLAST_NK4$V1 %in% as.character(dupdf_H$V1),]
noMultimatch_NK4$V11_V2 <- paste0(noMultimatch_NK4$V11,"|", noMultimatch_NK4$V2)
noMultimatch_NK4 <- noMultimatch_NK4[order(noMultimatch_NK4$V11),]
NK4_uniq <- as.data.frame(unique(noMultimatch_NK4$V1))
colnames(NK4_uniq) <- c("V1")
NK4_uniq <- join(NK4_uniq, noMultimatch_NK4, by = "V1", match = "first")

## Add marker location on the contigs
NK4_uniq$ID <- paste0(NK4_uniq$V1,"_",NK4_uniq$V2)
blast_NK4$ID <- paste0(blast_NK4$V1,"_",blast_NK4$V2)
blast_NK4 <- blast_NK4[order(blast_NK4$V9),]
NK4_uniq <- join(NK4_uniq, blast_NK4[,c("V9","ID")], by = "ID", match = "first")
blast_NK4 <- blast_NK4[order(blast_NK4$V10, decreasing = T),]
NK4_uniq <- join(NK4_uniq, blast_NK4[,c("V10","ID")], by = "ID", match = "first")
names(NK4_uniq)
colnames(NK4_uniq)[7:8] <- c("start", "stop")
colnames(NK4_uniq)[2]<- "NK4_Marker"

# Add LG and cM to matched markers
NK4_new <- merge(NK4_uniq, NK4[,c(2,5,10)],by="NK4_Marker")
NK4_new <- NK4_new[,c(1,3,9,10)]
colnames(NK4_new)[2:4] <- c("contig","NK4_Group","NK4_cM")


#-------------------oak map blast filtering-------------------------------------------
blast_Oak <- read.table("./Blast_Mummer_postChimeraRemoval/Oak_Cm.txt", sep = '\t')

## remove multiple matches, if a marker mapped to different contigs with same evalue, we discard it.
filterBLAST_Oak <- aggregate(V11~ V1+V2, data= blast_Oak, FUN=min)
filterBLAST_Oak <- filterBLAST_Oak[order(filterBLAST_Oak$V1),]
filterBLAST_Oak$V1_V11 <- paste0(filterBLAST_Oak$V1,".",filterBLAST_Oak$V11)
rows_dup <- which(filterBLAST_Oak$V1_V11 %in% unique(filterBLAST_Oak$V1_V11[duplicated(filterBLAST_Oak$V1_V11)]),)
dupdf <- filterBLAST_Oak[rows_dup,]
filterBLAST_Oak$V1 <- as.factor(filterBLAST_Oak$V1)
noMultimatch_Oak <- filterBLAST_Oak[!filterBLAST_Oak$V1 %in% as.character(dupdf$V1),]
noMultimatch_Oak <- noMultimatch_Oak[order(noMultimatch_Oak$V11),]
Oak_uniq <- as.data.frame(unique(noMultimatch_Oak$V1))
colnames(Oak_uniq) <- c("V1")
Oak_uniq <- join(Oak_uniq, noMultimatch_Oak, by = "V1", match = "first")

## Add marker location on the contigs
Oak_uniq$ID <- paste0(Oak_uniq$V1,"_",Oak_uniq$V2)
blast_Oak$ID <- paste0(blast_Oak$V1,"_",blast_Oak$V2)
blast_Oak <- blast_Oak[order(blast_Oak$V9),]
Oak_uniq <- join(Oak_uniq, blast_Oak[,c("V9","ID")], by = "ID", match = "first")
blast_Oak <- blast_Oak[order(blast_Oak$V10, decreasing = T),]
Oak_uniq <- join(Oak_uniq, blast_Oak[,c("V10","ID")], by = "ID", match = "first")
names(Oak_uniq)
colnames(Oak_uniq)[6:7] <- c("start", "stop")

############### Combine maps#####################################
## Merge maps
#1. combine mummer and blast to peach
peach2chestnut <- uniq_contigs[,c(1,7,8,9)]
names(peach2chestnut)
colnames(peach2chestnut) <- c("Location", "Tags", "Peach", "V1")
Pp2_CM_blast <- merge(peach2chestnut,peach_uniq, by = "V1", all = T)
names(Pp2_CM_blast)
Pp2_CM_blast <- Pp2_CM_blast[,-c(3,5,7,8)]
names(Pp2_CM_blast)
colnames(Pp2_CM_blast)[1] <- c("contig")

# Add orientations for Peach to chestnut using mummer and blast
#1. Mummer and BLAST directions are saved in "MummDir.csv", first three columns are mummer result, the rest are blast results.
#Pp2_CM_dir <- read.csv("./Blast_Mummer_postChimeraRemoval/MummDir.csv",header = T)
#head(Pp2_CM_dir)
#head(uniq_contigs)
#uniq_contigs_dir <- merge(uniq_contigs, Pp2_CM_dir[,c(1,3)], by.x = "Tag_Cm", by.y = "Cm", all.x = T)
#colnames(uniq_contigs_dir)[1] <- "contig"
#Pp2_CM_blast_dir <-join(Pp2_CM_blast, uniq_contigs_dir[,c(1,10)], by = "contig", type = "left", match = "first") 
#colnames(Pp2_CM_dir)[4] <- "contig"
#Pp2_CM_blast_dir <- join(Pp2_CM_blast_dir, Pp2_CM_dir[,c(4,8)], by = "contig", type = "left", match = "first")

##---
# need to recreate Pp2_CM_blast_dir
##---
###          contig            Location Peach   V2    start     stop   V13 Direction
### 1 contig0000003                <NA>  <NA> Pp04  2618049 25237118 minus        NA
### 2 contig0000004  9303536  9303758    Pp01 Pp01  2794855 34482605 minus        -1
#2. Add Fan's map
colnames(Fan_uniq)[3] <- c("contig")
Pp2_CM_blast_dir_Fan <- merge(Pp2_CM_blast_dir,Fan_uniq[,c(2,3)], by = "contig", all = T)
names(Pp2_CM_blast_dir_Fan)
colnames(Pp2_CM_blast_dir_Fan)[9] <- paste0("Markers")
names(Pp2_CM_blast_dir_Fan)
Pp2_CM_Fan <- join(Pp2_CM_blast_dir_Fan, FanMap, by = "Markers", match = "first")
names(Pp2_CM_Fan)[names(Pp2_CM_Fan) %in% c("Markers","cM","LG")]<-c("Fan_marker","Fan_cM","Fan_LG")

#3. Add Oak map
names(Oak_uniq)
colnames(Oak_uniq)[3] <- c("contig")
Pp2_CM_Fan_Oak <- merge(Pp2_CM_Fan, Oak_uniq[,c(2,3)], by = "contig", all = T)
names(Pp2_CM_Fan_Oak)
colnames(Pp2_CM_Fan_Oak)[12] <- c("Marker")
Pp2_CM_Fan_Oak <- join(Pp2_CM_Fan_Oak, Oak, by = "Marker", match="first")
names(Pp2_CM_Fan_Oak)[names(Pp2_CM_Fan_Oak) %in% c("Marker","LG","Position")]<-c("Oak_marker","Oak_LG","Oak_cM")

#4. Add HB2, JB1 and NK4
newMap <- merge(Pp2_CM_Fan_Oak,HB2_new, by="contig", all = T)
newMap_JB1 <- merge(newMap, JB1_new, by= "contig", all = T)
newMap_all <- merge(newMap_JB1, NK4_new, by = "contig", all = T)

# Rename some columns
names(newMap_all)
names(newMap_all)[names(newMap_all) %in% c("Peach","V2","V13","Direction")]<-c("Mumm_peach","Blst_peach","Blst_dir","Mumm_dir")

# Change direction labels to be consistent in mummer and blast.
newMap_all$Blst_dir <- gsub("minus","-1",newMap_all$Blst_dir)
newMap_all$Blst_dir <- gsub("plus","1", newMap_all$Blst_dir)

#5 Add Kubisiak map
colnames(Kub)[1:3] <- c("Fan_marker","Kub_LC","Kub_cM")
newMap_all_0223 <- join(newMap_all, Kub, by = "Fan_marker", match = "first")

# split the location column
library(tidyr)
newMap_all_0223$Location <- lapply(newMap_all_0223$Location, trimws)
newMap_all_0223 <- separate(data=newMap_all_0223, col = Location, into = c("Mummer_peach_start", "Mummer_peach_end"), sep = "\\s+")

# fix column names and order
colnames(newMap_all_0223)[c(1,2,6)] <- c("Kub_Fan_marker", "Chestnut_contig", "Mummer_peach_scaffold")
colnames(newMap_all_0223)[c(7,8,9)]<- c("Blast_peach_scaffold", "Blast_peach_start", "Blast_peach_end")
colnames(newMap_all_0223)[c(10,11,17)]<- c("Blast_direction", "Mummer_direction", "HB2_LG")
colnames(newMap_all_0223)[c(20,23,25)]<- c("JB1_LG", "NK4_LG", "Kub_LG")

newMap_all_0223 <- newMap_all_0223[c("Chestnut_contig", "Mummer_peach_scaffold", "Mummer_peach_start", "Mummer_peach_end", "Mummer_direction", "Blast_peach_scaffold", "Blast_peach_start", "Blast_peach_end", "Blast_direction", "Kub_Fan_marker", "Fan_LG", "Fan_cM", "Kub_LG", "Kub_cM", "HB2_Marker", "HB2_LG", "HB2_cM", "JB1_Marker", "JB1_LG", "JB1_cM", "NK4_Marker", "NK4_LG", "NK4_cM", "Oak_marker", "Oak_LG", "Oak_cM")]

write.csv(newMap_all_0223,"AllMaps_20180223.csv")

#write.csv(newMap_all,"AllMaps_20180212.csv")
### Need to split the mummer reference location and change to order of columns in excel sheet. 

#------------count duplications in Fan edited map----------------------------
# Fan's 2nd edit, remove some duplications
Fan_edit2 <- read.csv("Chinese_only_integrated_edit2.csv", header = T)
colnames(Fan_edit2)[3] <- "LG"
Fan_edit2 <- Fan_edit2[-c(grep("group",Fan_edit2$Markers)),]
Fan_edit2$Markers <- gsub("_\\d+","\\_v2", Fan_edit2$Markers)
Fan_edit2$Markers <- gsub("v2c","\\_contig", Fan_edit2$Markers)
duplicates <- Fan_edit2[c(which(Fan_edit2$Markers %in% unique(Fan_edit2$Markers[duplicated(Fan_edit2$Markers)]),)),]

