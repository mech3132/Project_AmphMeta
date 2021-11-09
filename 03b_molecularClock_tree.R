#!bin/bash RScript

library(ape)
library(tidyverse)

dir.create("03a_molecularClock_tree/Routput")
dir.create("03a_molecularClock_tree/Routput/pathd8_nodes")

# 
tree <- read.tree("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/tree_99_ASV50/tree.nwk")
# treeUM <- read.tree("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/tree_99_ASV50/tree_UM.nwk")
taxa <- read.delim("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/taxonomy_vsearch_80/taxonomy.tsv")
silva <- read.delim("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/databases/silva-138-99-tax-515-806/taxonomy.tsv")
otu <- read.delim("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt/feature-table.txt", skip=1, row.names=1)
# Get taxonomies to filter
taxa <- taxa %>% 
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species")
           , sep="; ", remove=FALSE, fill="right") %>%
  mutate(Phylum = ifelse(is.na(Phylum), Domain, Phylum)
         , Class = ifelse(is.na(Class), Phylum, Class)
         , Order = ifelse(is.na(Order), Class, Order)
         , Family = ifelse(is.na(Family), Order, Family)
         , Genus = ifelse(is.na(Genus), Family, Genus)
         , Species = ifelse(is.na(Phylum), Genus, Species)) 
silva <- silva %>% filter(Feature.ID %in% tree$tip.label) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species")
           , sep="; ", remove=FALSE, fill="right") %>%
  mutate(Phylum = ifelse(is.na(Phylum), Domain, Phylum)
         , Class = ifelse(is.na(Class), Phylum, Class)
         , Order = ifelse(is.na(Order), Class, Order)
         , Family = ifelse(is.na(Family), Order, Family)
         , Genus = ifelse(is.na(Genus), Family, Genus)
         , Species = ifelse(is.na(Phylum), Genus, Species)) 
# Get Archaea
archaea_id <- silva %>% filter(Domain == "d__Archaea") %>%
  group_by(Domain) %>% mutate(rank=rank(Feature.ID)) %>% filter(rank==1) %>% pull(Feature.ID)

# Get cyanobacteria
cyano_ids <- silva %>% filter(Phylum=="p__Cyanobacteria") %>% group_by(Order) %>%
  mutate(rank = rank(Feature.ID)) %>% filter(rank==1) %>% pull(Feature.ID)

# Get ricketts
ricket_ids <- silva %>% filter(Family=="f__Rickettsiales") %>% group_by(Genus) %>%
  mutate(rank = rank(Feature.ID)) %>% filter(rank<=3) %>% pull(Feature.ID)

# Get Chlorobium
chloro_ids <- silva %>% filter(Class=="c__Chlorobia") %>% group_by(Genus) %>%
  mutate(rank = rank(Feature.ID)) %>% filter(rank==1) %>% pull(Feature.ID)

# Get Chromat
chromat_ids <- silva %>% filter(Order=="o__Chromatiales") %>% group_by(Genus) %>%
  mutate(rank = rank(Feature.ID)) %>% filter(rank==1) %>% pull(Feature.ID)

# Get all tips
allMolecClockIDs <- c(cyano_ids, ricket_ids, chloro_ids, chromat_ids, archaea_id)
## Add in some extras
# topOTUs <- names(sort(rowSums(otu))[1:100])
tree_molecOnly <- keep.tip(tree, c(allMolecClockIDs))
# tree_molecOnlyUM <- keep.tip(treeUM, c(allMolecClockIDs,topOTUs))

tree_meta <- data.frame(ASVID=tree_molecOnly$tip.label) %>% 
  mutate(group=ifelse(ASVID %in% cyano_ids, "Cyanobacteria", ifelse(ASVID %in% ricket_ids, "Rickettsiales", ifelse(ASVID %in% chromat_ids, "Chromatiales", ifelse(ASVID %in% chloro_ids, "Chlorobi", "Archaea"))))) %>%
  mutate(color= ifelse(group=="Cyanobacteria", "turquoise", ifelse(group=="Rickettsiales", "brown", ifelse(group== "Chromatiales", "red", ifelse(group=="Chlorobi", "green", "black")))))
# tree_molecOnly$tip.label <- paste0(silva[match(tree_molecOnly$tip.label, silva$Feature.ID), "Class"], tree_molecOnly$tip.label)
png(filename="03a_molecularClock_tree/Routput/molecular_clock_sequence_tree.png", height = 700, width=1000)
plot(tree_molecOnly, tip.color = tree_meta$color)
dev.off()
# png(filename="03a_molecularClock_tree/Routput/molecular_clock_sequence_treeUM.png", height = 700, width=1000)
# plot(tree_molecOnlyUM, tip.color = tree_meta$color)
# dev.off()

### Get filtered tree
allTips <- c(rownames(otu), allMolecClockIDs)
tree_filt <- keep.tip(tree, allTips)

### Go through molecular clock dating
cyano_table <- c()
for ( i1 in 1:(length(cyano_ids)-1) ) {
  for (i2 in 2:length(cyano_ids)) {
    cyano_table <- c(cyano_table, paste0(c("mrca: ", cyano_ids[i1], ", ",cyano_ids[i2], ", minage=0.0025;"), collapse=""))
  }
}
ricket_table <- c()
for ( i1 in 1:(length(ricket_ids)-1) ) {
  for (i2 in 2:length(ricket_ids)) {
    ricket_table <- c(ricket_table, paste0(c("mrca: ", ricket_ids[i1], ", ",ricket_ids[i2], ", minage=1.6;"), collapse=""))
  }
}

chromat_table <- c()
for ( i1 in 1:(length(chromat_ids)-1) ) {
  for (i2 in 2:length(chromat_ids)) {
    chromat_table <- c(chromat_table, paste0(c("mrca: ", chromat_ids[i1], ", ",chromat_ids[i2], ", minage=1.64;"), collapse=""))
  }
}

chloro_table <- c()
for ( i1 in 1:(length(chloro_ids)-1) ) {
  for (i2 in 2:length(chloro_ids)) {
    chloro_table <- c(chloro_table, paste0(c("mrca: ", chloro_ids[i1], ", ",chloro_ids[i2], ", minage=1.64;"), collapse=""))
  }
}
# Get divergence of archaea to everything else
archaea_table <- c()
for ( bact in allMolecClockIDs ) {
  archaea_table <- c(archaea_table, paste0(c("mrca: ", bact, ", ",archaea_id, ", fixage=3.8;"), collapse=""))
}


seqString <- paste0("Sequence length = 150;")

write.tree(tree_filt, file="03a_molecularClock_tree/Routput/filtered_tree.nwk")
write.table(seqString, "03a_molecularClock_tree/Routput/seqString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)
write.table(data.frame(cyano_table), "03a_molecularClock_tree/Routput/pathd8_nodes/cyanoString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)
write.table(data.frame(ricket_table), "03a_molecularClock_tree/Routput/pathd8_nodes/ricketString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)
write.table(data.frame(chloro_table), "03a_molecularClock_tree/Routput/pathd8_nodes/chloroString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)
write.table(data.frame(chromat_table), "03a_molecularClock_tree/Routput/pathd8_nodes/chromatString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)
write.table(data.frame(archaea_table), "03a_molecularClock_tree/Routput/pathd8_nodes/archaeaString.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)

# TEST
write.tree(tree_molecOnly, file="03a_molecularClock_tree/Routput/tree_Silva.nwk")
# write.tree(tree_molecOnlyUM, file="03a_molecularClock_tree/Routput/tree_Silva_UM.nwk")


## COMPARE
# calTree <- read.tree("03a_molecularClock_tree/final_tree_99_ASV50.nwk")
# 
# topOTUs <- names(sort(rowSums(otu))[1:100])
# calTree_top <- keep.tip(calTree, topOTUs)
# name_taxa <- make.unique(taxa[match(calTree_top$tip.label, taxa$Feature.ID),"Class"])
# calTree_top$tip.label <- name_taxa
# plot(calTree_top )


