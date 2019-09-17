## get stats for the psudochromasome table

# read tables
new_ps <- read.csv("Pseudochromosomes_cdn_bac_v11.csv", header = T)
# read the contig length table
contig_len <- read.table("Castanea_mollissima_scaffolds_v3.2.lens.tsv", header = F)
names(contig_len) <- c("cc_contig", "length")
contig_len$cc_contig <- gsub("lcl\\|","", contig_len$cc_contig)
# merge the lengths to psudochro table
new_ps <- join(new_ps, contig_len, by = "cc_contig", type="left")

# get gene numbers for each v3.2 contigs
gff <- read.table("Castanea_mollissima_scaffolds_v3.2_ALLmRNA.gff")
dedup_gff <- gff[grep("t1", gff$V9),]
contig_gene_number <- data.frame(table(dedup_gff$V1))
names(contig_gene_number) <- c("cc_contig", "gene_number")
contig_gene_number$cc_contig <- gsub("lcl_","", contig_gene_number$cc_contig)

# merge gene numbers to the psudochro table
new_ps <- join(new_ps, contig_gene_number, by = "cc_contig", type = "left")
# seperate multimap contigs and unique map contigs
rows_dup <- which(new_ps$cc_contig %in% unique(new_ps$cc_contig[duplicated(new_ps$cc_contig)]),)
multimap <- new_ps[rows_dup,]
length(unique(multimap$cc_contig))
sum(unique(multimap$length))
uniq <- new_ps#[-rows_dup,]

sum(uniq$length,na.rm=T)
length(uniq$cc_contig)
sum(multimap$gene_number[unique(multimap$cc_contig)], na.rm = T)
sum(uniq$length, na.rm = T)
sum(uniq$gene_number, na.rm = T)
multimap_uniq <- new_ps[unique(multimap$cc_contig),]
length(which(is.na(multimap_uniq$gene_number)))
length(which(is.na(uniq$gene_number)))

# count bac_ctg numbers
bac <- new_ps[grep("bac",new_ps$cc_contig),]
length(unique(bac$cc_contig))

# count v4.1 gene numbers
names(contig_gene_number) <- c("Chestnut_contig", "gene_number")
CmVersion4 <- join(CmVersion4, contig_gene_number, by = "Chestnut_contig", type = "left")

# save contig gene number into a csv file
write.csv(contig_gene_number, "Cm_contigs_gene.numbers.csv") 
sum(contig_gene_number$gene_number, rm.ma=T)
