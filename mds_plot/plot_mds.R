setwd("~/Desktop/Jiali/UTK/chestnut/chestnut_genome/mds_plot/")
library(ggplot2)
## MDS plot using all.fitered.mds
d <- read.table("plink.mds", h=T)
d$genotype <-  gsub(":H.*","",d$FID)
d$genotype[1:76] <- "BC3F3"
d <- d[-(191:192),]
d$species = factor(c(rep("BC3F3", 76), rep("BC1 Clapper", 16), rep("C.dentata", 82), rep("Vanuxem",16)))

ggplot(d, aes(C1, C2, color=genotype, shape=species)) +
  geom_point(size=3) +
  xlab("MDS 1") +
  ylab("MDS 2") + 
  coord_fixed() +
  scale_color_brewer(palette="Set3")+
  theme_classic(base_size = 14) +
  guides(size = guide_legend(order=2),
         shape = guide_legend(order=1))


ggsave("MDS_all_version3.2.png", height = 8, width = 8)  

## mds plot using all.taxamerge.filtered.mds
taxa <- read.table("all.taxamerged.filtered.mds", h=T, stringsAsFactors = F)
taxa$IID[1:38] <- "B3F3"
names(taxa)[2] <- 'genotype'
taxa$species = factor(c(rep("B3F3", 38), "C.dentata", "Clapper", rep("C.dentata", 8), "Vanuxem"))

ggplot(taxa[-c(43,44),], aes(C1, C2, shape=species, color=genotype)) +
  geom_point(size=3) +
  xlab("MDS component 1") +
  ylab("MDS component 2") + 
  coord_fixed() +
  scale_color_brewer(palette="Set3")+
  theme_classic(base_size = 14)

ggsave("MDS_taxa_merged_excluded.png", height = 8, width = 8)
