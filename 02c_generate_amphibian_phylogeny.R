#!bin/bash Rscript
library(ape)
library(tidyverse)

amphtree <- read.tree("00_amphibian_phylogeny/amph_shl_new_Consensus_7238.tre")
# amphCSV <- read.csv("00_amphibian_phylogeny/amph_shl_new_Classification.csv")
# meta <- read.delim("01_parsed_meta/full_meta.txt")
meta <- read.delim("02b_bioclim_ecology/meta_complete.txt")
amphfam <- read.csv("00_amphibian_phylogeny/amph_shl_new_Classification.csv")
# amphecoRegion <- read.csv("00_amphibian_phylogeny/Global_Ecoregion_Assignments.csv")


dir.create("02c_generate_amphibian_phylogeny")

#####
# Get ecoregion legend
# amphecoRegion_legend <- amphecoRegion %>% select(Areas) %>% 
#   mutate(Areas=ifelse(Areas=="",NA,Areas)) %>%
#   drop_na() %>%
#   separate(Areas, into=c("Number","EcoRegion"), sep ="-") %>%
#   mutate(Number=as.numeric(Number))

# Match Scientific name with tree tips
amphfam_adj <- amphfam %>% mutate(HostPhylo = gsub(" ","_",Scientific.Name)) %>%
  rename(AmphTaxon=Taxon, AmphFamily=Family, AmphSubFamily=Subfamily, AmphClade=Clade) %>%
  select(HostPhylo, AmphTaxon, AmphFamily, AmphSubFamily, AmphClade)
# amphecoRegion_adj <- amphecoRegion %>% rename(HostPhylo=Species, Number=Region) %>%
#   mutate(Number=as.numeric(Number)) %>% 
#   left_join(amphecoRegion_legend) %>% select(HostPhylo, EcoRegion)
# amphecoRegion_adj %>% mutate(duplicated=duplicated(HostPhylo)) %>% filter(duplicated) %>%
#   select(HostPhylo) %>% left_join(amphecoRegion_adj)

meta_adj <- meta %>% rowwise() %>%mutate(HostPhylo = gsub(" ","_", Host), HostOriginal=Host)  %>%
  mutate(HostPhylo = ifelse(HostPhylo == "Boana_faber", "Hypsiboas_faber",
                            ifelse(HostPhylo == "Boana_albopunctata", "Hypsiboas_albopunctatus",
                                   ifelse(HostPhylo == "Boana_bischoffi", "Hypsiboas_bischoffi",
                                          ifelse(HostPhylo == "Boana_semilineata", "Hypsiboas_semilineatus",
                                                 ifelse(HostPhylo == "Boana_albomarginata", "Hypsiboas_albomarginatus",
                                                        ifelse(HostPhylo == "Boana_polytaenia", "Hypsiboas_polytaenius",
                                                               ifelse(HostPhylo == "Boana_atlantica", "Hypsiboas_atlanticus", 
                                                                      ifelse(HostPhylo == "Dendrosophus_ebraccatus", "Dendropsophus_ebraccatus", 
                                                                                    ifelse(HostPhylo == "Glandirana_rugosa", "Rugosa_rugosa", HostPhylo
                                   )))))))))) %>%
  mutate(Host = gsub(" ", "_", Host)) %>% ungroup() %>%
  left_join(amphfam_adj)
# nrow(meta)
# nrow(meta_adj)


allSpecies <- unique(meta_adj$HostPhylo)
allSpecies[which(!allSpecies %in% amphtree$tip.label)]

amphtree_filt <- keep.tip(amphtree, tip=allSpecies)
# amph_dist <- cophenetic(amphtree_filt)
amphtree_filt$tip.label <- meta_adj$Host[match(amphtree_filt$tip.label, meta_adj$HostPhylo)]

#### Save amphibian
write.tree(amphtree_filt, file="02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
# write.table(amph_dist, )

## Adjust meta headers for . vs - conversion

write.table(meta_adj, file = "02c_generate_amphibian_phylogeny/metadata_hostcorrected.txt", quote=FALSE, row.names=FALSE, sep="\t")

png("02c_generate_amphibian_phylogeny/amphibian_phylogeny.png", height=1250, width=750)
plot(amphtree_filt)
dev.off()
#### With ecotype

HabClass <- meta_adj[match(amphtree_filt$tip.label, meta_adj$Host), "HabitatClass"]
colorMatch <- HabClass %>% distinct() %>% mutate(colour=ifelse(HabitatClass=="Aquatic","blue"
                                                 ,ifelse(HabitatClass=="Semi-aquatic","turquoise"
                                                         , ifelse(HabitatClass=="Terrestrial","brown", 
                                                                  ifelse(HabitatClass=="Arboreal","green", "black")))))
colorTips <- HabClass %>% left_join(colorMatch) %>% pull(colour)
png("02c_generate_amphibian_phylogeny/amphibian_phylogeny_withEco.png", height=1250, width=750)
plot(amphtree_filt, tip.color = colorTips)
dev.off()
