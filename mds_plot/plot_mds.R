setwd("~/Desktop/Jiali/UTK/chestnut/chestnut_genome/mds_plot/")
library(ggplot2)
## MDS plot using all.fitered.mds
d <- read.table("all.filtered.mds", h=T)
d$genotype <-  gsub(":H.*","",d$FID)
d$genotype[1:74] <- "B3F3"
d$species = factor(c(rep("B3F3", 74), rep("Clapper", 16), rep("C.dentata", 79), rep("Vanuxem",16)))

plot(d$C1, d$C2, col=as.integer(d$species), pch=19, xlab="MDS component 1", ylab="MDS component 2", main = "MDS")
legend("topright", c("B3F3", "Clapper", "C.dentata","Vanuxem"), pch=19, col=c(1,3,2,4))

ggplot(d, aes(C1, C2, shape=species, color=genotype)) +
  geom_point(size=3) +
  xlab("MDS component 1") +
  ylab("MDS component 2") + 
  coord_fixed() +
  scale_color_brewer(palette="Set3")+
  theme_classic(base_size = 14)

ggsave("MDS_all.png", height = 8, width = 8)  

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
