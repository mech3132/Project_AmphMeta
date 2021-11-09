#!bin/bash
library(ape)
library(maptools)
library(maps)
library(tidyverse)
dir.create("04e_partition_data_for_downstream")

#### Load data #####
# 
# otu_rare <- read.delim("04b_data_filt_and_combine/downstream/otu_rare.txt", row.names=1, sep = "\t")
# otu_cut <- read.delim("04b_data_filt_and_combine/downstream/otu_cut.txt", row.names = 1, sep = "\t")

meta_cut <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt", sep = "\t")
## Filter Hyla Japonica this one is weird
meta_cut <- meta_cut %>% filter(Host !="Hyla_japonica")
otu_cut <- read.delim("04b_data_filt_and_combine/downstream/otu_cut.txt")
# 
# bcdm <- read.delim("04b_data_filt_and_combine/downstream/dm_bray.txt", row.names = 1, sep="\t")
# wudm <- read.delim("04b_data_filt_and_combine/downstream/dm_wu.txt", row.names = 1, sep="\t")
# uwudm <- read.delim("04b_data_filt_and_combine/downstream/dm_uwu.txt", row.names = 1, sep="\t")

amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")


#### Abundant hosts only ####
# highPrevTaxa <- meta_cut %>% 
#   group_by(Host, Location, CollectionDate) %>%
#   summarize(totalUnique=n()) %>% ungroup() %>%
#   filter(totalUnique>5, !is.na(CollectionDate)) %>%
  # group_by(Host) %>% summarize(nUnique=n()) %>%
  # filter(nUnique>=3) %>% ungroup() %>%
  # pull(Host)
toKeepAbundant <- meta_cut %>% #filter(Host %in% highPrevTaxa, !is.na(CollectionDate)) %>% select(Host, Location, CollectionDate) %>%
  filter(!is.na(CollectionDate), !is.na(Host), !is.na(Latitude), !is.na(Longitude), !is.na(HabitatClass), !is.na(Location)) %>%
  group_by(Host,Location,CollectionDate) %>% summarize(perBout=n()) %>% ungroup() %>%
  arrange(perBout) %>% filter(perBout>5) %>%
  arrange(Host) %>%
  select(Host, Location, CollectionDate) %>% distinct()
# View(toKeepAbundant)

meta_cut %>% filter(Host %in% toKeepAbundant$Host, Location %in% toKeepAbundant$Location, CollectionDate %in% toKeepAbundant$CollectionDate) %>%
  nrow()

meta_cut %>% filter(Host %in% toKeepAbundant$Host, Location %in% toKeepAbundant$Location, CollectionDate %in% toKeepAbundant$CollectionDate) %>%
  select(Host, Location, Latitude, Longitude, HabitatClass, Month) %>%
  drop_na() %>% nrow()

meta_final <- meta_cut %>% filter(Host %in% toKeepAbundant$Host, Location %in% toKeepAbundant$Location, CollectionDate %in% toKeepAbundant$CollectionDate) %>%
  select(sampleid, Study, ENASTUDY, Host, Location, Latitude, Longitude, HabitatClass, Month, CollectionDate, starts_with("Amph"), starts_with("biom_")) %>%
  drop_na() 
nrow(meta_final)

write.table(toKeepAbundant, file="04e_partition_data_for_downstream/abundant_amphibians.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(meta_final, file="04e_partition_data_for_downstream/final_meta.txt", row.names = FALSE, quote = FALSE, sep="\t")

# Filter OTU table
otu_filt <- otu_cut[,meta_final$sampleid]
otu_filt <- otu_filt[rowSums(otu_filt)>0,]
nrow(otu_filt)
sum(otu_filt)

sink("04e_partition_data_for_downstream/abundant_amphibian_notes.txt")
print("unique Host/Loc/Date greater than 5")
print("Number of samples")
meta_final %>% nrow()
print("Number of species")
meta_final %>% select(Host) %>% pull() %>% unique() %>% length()
print("Number of studies" )
meta_final %>% select(Study) %>% pull() %>% unique() %>% length()
print("Number of OTUs")
nrow(otu_filt)
print("Number of reads")
sum(otu_filt)

sink()



# 
# # Quick plot
# meta_abundant <- meta_cut %>% 
#   filter(Host %in% toKeepAbundant$Host,Location %in% toKeepAbundant$Location, CollectionDate %in% toKeepAbundant$CollectionDate) 
# gg_abundHost <- meta_abundant %>%
#   group_by(Host, Location) %>% summarize(nSamples=n()) %>% ungroup() %>%
#   ggplot() + geom_point(aes(x=Host, y=Location, col=log10(nSamples)), cex=4) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +xlab("Host species")
# gg_abundHost
# # ggsave(filename = "04b_data_filt_and_combine/abundant_host_samples.png", height=5, width=8, gg_abundHost)
# keepTipCol <- rep("black", length(amphtree$tip.label))
# keepTipCol[match(meta_abundant$Host, amphtree$tip.label)] <- "red"
# plot(amphtree, tip.col=keepTipCol)
# plot(keep.tip(amphtree, unique(meta_abundant$Host)))
# length(unique(meta_abundant$Host))
# length(unique(meta_abundant$Study))
# 
# # png("04d_host_and_enviro_matrices/map_samples_sizeSamp_colRich.png", height=600, width=900)
# maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
# points(meta_abundant$Longitude,meta_abundant$Latitude,bg="red", pch=21)
# # dev.off()
# 
# 
# #### TEMPORARY SMALL DATASET FOR WRITING CODE
# 
# toKeepAbundant2 <- meta_cut %>% #filter(Host %in% highPrevTaxa, !is.na(CollectionDate)) %>% select(Host, Location, CollectionDate) %>%
#   filter(!is.na(CollectionDate)) %>%
#   filter(!is.na(Latitude), !is.na(Longitude)) %>%
#   group_by(Host,Location,CollectionDate) %>% summarize(perBout=n()) %>% ungroup() %>%
#   arrange(perBout) %>% filter(perBout>20) %>%
#   arrange(Host) %>%
#   select(Host, Location, CollectionDate) %>% distinct() 
# # View(toKeepAbundant)
# 
# nrow(meta_cut)
# meta_cut %>% filter(Host %in% toKeepAbundant2$Host, Location %in% toKeepAbundant2$Location, CollectionDate %in% toKeepAbundant2$CollectionDate) %>%
#   nrow()
# 
# write.table(toKeepAbundant2, file="04e_partition_data_for_downstream/abundant_amphibians_small.txt", row.names = FALSE, quote = FALSE, sep="\t")
# 
# 
# # Quick plot
# meta_abundant2 <- meta_cut %>% 
#   filter(Host %in% toKeepAbundant2$Host,Location %in% toKeepAbundant2$Location, CollectionDate %in% toKeepAbundant2$CollectionDate) 
# gg_abundHost <- meta_abundant2 %>%
#   group_by(Host, Location) %>% summarize(nSamples=n()) %>% ungroup() %>%
#   ggplot() + geom_point(aes(x=Host, y=Location, col=log10(nSamples)), cex=4) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +xlab("Host species")
# gg_abundHost
# # ggsave(filename = "04b_data_filt_and_combine/abundant_host_samples.png", height=5, width=8, gg_abundHost)
# keepTipCol <- rep("black", length(amphtree$tip.label))
# keepTipCol[match(meta_abundant2$Host, amphtree$tip.label)] <- "red"
# plot(amphtree, tip.col=keepTipCol)
# plot(keep.tip(amphtree, unique(meta_abundant2$Host)))
# length(unique(meta_abundant2$Host))
# length(unique(meta_abundant2$Study))
# 
# # png("04d_host_and_enviro_matrices/map_samples_sizeSamp_colRich.png", height=600, width=900)
# maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
# points(meta_abundant2$Longitude,meta_abundant2$Latitude,bg="red", pch=21)
# # dev.off()
# 
# ###### 
# # Species that span multiple studies
# multipleStudyLocList <- meta_cut %>%
#   group_by(Host) %>%
#   summarise(nStudy=length(unique(Study)), nLoc=length(unique(Location))) %>%
#   filter(nStudy>1, nLoc>1)
# write.table(multipleStudyLocList, file="04e_partition_data_for_downstream/multipleStudyLocList.txt", row.names = FALSE, quote = FALSE, sep="\t")
# 
# # SEe amphibians
# plot(keep.tip(amphtree, multipleStudyLocList$Host))
# 
# ### Studies with multiple locations and species
# # Look for effects of location and species
# meta_cut %>%
#   unite(Host, Location, col="SpLoc", remove=FALSE) %>%
#   group_by(Study) %>%
#   summarize(nLoc = length(unique(Location)), nSp = length(unique(Host))
#             , nLocxSp=length(unique(SpLoc)), ntotal=n()) %>%
#   filter(nLoc>3, nSp>3)
# 
# ## Species at multiple locations and studies
# Species_multipleLocStudy <- meta_cut %>%
#   group_by(Host) %>%
#   summarize(nLoc = length(unique(Location)), nStudy=length(unique(Study))) %>%
#   filter(nStudy>=3, nLoc>nStudy) %>% pull(Host)
# 

