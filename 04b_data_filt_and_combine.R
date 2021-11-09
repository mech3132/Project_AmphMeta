#!bin/bash

library(MASS) #NMDS
library(lubridate)
library(ape) # to get the Branch * sites matrices
library(tidyverse)

#### Loading ####
meta <- read.delim("02c_generate_amphibian_phylogeny/metadata_hostcorrected.txt", sep="\t")
otu_raw <- read.delim("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt/feature-table.txt", skip=1, row.names=1)
otu <- read.delim("04a_diversity/filtered_tables/table_cut2500/feature-table.txt", skip=1, row.names = 1)
otu_rare <- read.delim("04a_diversity/filtered_tables/table_rare2500/feature-table.txt", skip=1, row.names = 1)
tree <- read.tree("04a_diversity/tree/tree_filt/tree.nwk")
tree_cal <- read.tree("04a_diversity/tree/calibrated_tree_filt/tree.nwk")
amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
taxa <- read.delim("02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/taxonomy_vsearch_80/taxonomy.tsv")
## Diversity load
obsotu <- read.delim("04a_diversity/exported/observed_otus.tsv") %>% rename(sampleid_original=X)
faithpd <- read.delim("04a_diversity/exported/faith_pd.tsv") %>% rename(sampleid_original=X)
bray <- read.delim("04a_diversity/exported/bray_curtis.tsv", row.names = 1, header = TRUE)
rownames(bray) <- colnames(bray)
# uwu <- read.delim("04a_diversity/exported/unweighted_unifrac.tsv", row.names = 1, header = TRUE)
# rownames(uwu) <- colnames(uwu)
wu <- read.delim("04a_diversity/exported/weighted_unifrac.tsv", row.names = 1, header = TRUE)
rownames(wu) <- colnames(wu)
jac <- read.delim("04a_diversity/exported/jaccard_matrix_cut.tsv", row.names = 1, header = TRUE)
rownames(jac) <- colnames(jac)
#### ADD IN JACCARD

dir.create("04b_data_filt_and_combine")
dir.create("04b_data_filt_and_combine/downstream")

sink(file="04b_data_filt_and_combine/LOG.txt")
print("SUMMARY OF EDA")
sink()
########## Basic EDA ##########

#### Fix metadata ####
source("code/date_mapping.R")
date_mapping <- meta %>% select(ENASTUDY, CollectionDate) %>% mutate(FixedDate=NA) %>% distinct() %>%
  rename(referenceDate=CollectionDate)

date_mapping_complete <- convertDates(date_mapping, Study = "ENASTUDY")
# date_mapping_complete %>% mutate(FixedDate = as.Date(FixedDate)) %>% View()
##R loads OTU names with . instead of -, so I'm changing all- in MF to .
meta_fix <- meta %>% rename(referenceDate=CollectionDate) %>%
  left_join(date_mapping_complete) %>% rename(CollectionDate=FixedDate) %>%
  mutate(sampleid_original=sampleid) %>% rowwise() %>%
  mutate(sampleid = gsub("-",".",sampleid, fixed = TRUE)) %>%
  mutate(sampleid = ifelse(ENASTUDY=="Prest", gsub("_",".", sampleid, fixed=TRUE), sampleid)) %>%
  ungroup() 
sink(file="04b_data_filt_and_combine/LOG.txt", append = TRUE)
print("\nMetadata has all - converted to .; and _ converted to . for Prest. This is due to R's loading quirks")
sink()


## What samples are in metaAND in otu table?
# NOTE: print how many samples NOT in otu
sink(file="04b_data_filt_and_combine/LOG.txt", append = TRUE)
print("\nNumber of samples that could not be loaded into OTU table:")
sum(!meta_fix$sampleid %in% colnames(otu_raw))
sink()
meta_fix[which(!meta_fix$sampleid %in% colnames(otu_raw)),] 
# Most of these appear to be randos; maybe sequencing failed. 

meta_filt <- meta_fix[which(meta_fix$sampleid %in% colnames(otu)), ]
otu <- otu[,which(colnames(otu) %in% meta_filt$sampleid)] # This should already be filtered
otu_raw <- otu_raw[,which(colnames(otu_raw) %in% meta_fix$sampleid)] 
otu_rare <- otu_rare[,which(colnames(otu_rare) %in% meta_fix$sampleid)] 

any(!rownames(otu) %in% taxa$Feature.ID) # This should be false
# prestids %in% meta_filt$sampleid
#### Taxonomy ####
# Cross-check ASV and taxonomies
taxa_expand <- taxa %>% filter(Feature.ID%in%rownames(otu)) %>%
  rename(FeatureID = Feature.ID) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species")
           , sep="; ", remove=FALSE, fill="right") %>%
  rowwise() %>%
  mutate(Phylum = ifelse(length(grep("uncultured|Uncultured|ncertae|unknown|Unknown|metagenome",Phylum))>0, NA, Phylum)
         , Class = ifelse(length(grep("uncultured|Unculturedncertae|unknown|Unknown|metagenome",Class))>0, NA, Class)
         , Order = ifelse(length(grep("uncultured|Unculturedncertae|unknown|Unknown|metagenome",Order))>0, NA, Order)
         , Family = ifelse(length(grep("uncultured|Unculturedncertae|unknown|Unknown|metagenome",Family))>0, NA, Family)
         , Genus = ifelse(length(grep("uncultured|Unculturedncertae|unknown|Unknown|metagenome",Genus))>0, NA, Genus)
         , Species = ifelse(length(grep("uncultured|Unculturedncertae|unknown|Unknown|metagenome",Species))>0, NA, Species)) %>% 
  ungroup() %>%
  mutate(Phylum = ifelse(is.na(Phylum), Domain, Phylum)
         , Class = ifelse(is.na(Class), Phylum, Class)
         , Order = ifelse(is.na(Order), Class, Order)
         , Family = ifelse(is.na(Family), Order, Family)
         , Genus = ifelse(is.na(Genus), Family, Genus)
         , Species = ifelse(is.na(Species), Genus, Species)) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Class = make.unique(Class)) %>%
  mutate(Order = make.unique(Order)) %>%
  mutate(Family = make.unique(Family)) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Species = make.unique(Species)) %>%
  rowwise() %>%
  mutate(Phylum = ifelse(length(grep("p_", Phylum))>0, gsub("[.].*$", "", Phylum), Phylum)
         ,Class = ifelse(length(grep("c_", Class))>0, gsub("[.].*$", "", Class), Class)
         ,Order = ifelse(length(grep("o_", Order))>0, gsub("[.].*$", "", Order), Order)
         ,Family = ifelse(length(grep("f_", Family))>0, gsub("[.].*$", "", Family), Family)
         ,Genus = ifelse(length(grep("g_", Genus))>0, gsub("[.].*$", "", Genus), Genus)
         ,Species = ifelse(length(grep("s_", Species))>0, gsub("[.].*$", "", Species), Species)
  ) %>% 
  mutate(ChlorOrMito = length(grep("chloroplast|mitochondr", Taxon, ignore.case = TRUE))>0
         , Archaea = length(grep("d__Archaea", Taxon, ignore.case = TRUE))>0
         , Unassigned = length(grep("Unassigned", Taxon, ignore.case = TRUE))>0) %>%
  ungroup()

# Check everything in OTU table is in taxa
any(!rownames(otu) %in% taxa_expand$FeatureID) # Should be false

# Check the trees
any(!tree$tip.label %in% rownames(otu)) # Should be false
any(!tree_cal$tip.label %in% rownames(otu)) # Should be false

#### Histograms ####
# Distribution of samples- to cut off
hist_reads <- data.frame(Reads=colSums(otu_raw)) %>% rownames_to_column(var="sampleid") %>%
  left_join(meta_fix) %>%
  ggplot() + geom_histogram(aes(x=(Reads), fill=ENASTUDY))+
  scale_x_log10() + ylab("Count") + xlab("Reads per sample")
hist_reads
ggsave(filename = "04b_data_filt_and_combine/histogram_reads.png",height=5, width=7, hist_reads)
# Number of samples lost with 2500 filtering
sink(file="04b_data_filt_and_combine/LOG.txt", append = TRUE)
print("Number of samples lost with 2500 cutoff:")
ncol(otu_raw) - ncol(otu)
sink()

#### Dm filtering ####
bray_filt <- bray[meta_filt$sampleid, meta_filt$sampleid]
wu_filt <- wu[meta_filt$sampleid, meta_filt$sampleid]
jac_filt <- jac[meta_filt$sampleid, meta_filt$sampleid]

#### Alpha combine ####
meta_filt <- meta_filt %>% left_join(obsotu) %>% left_join(faithpd)

### Beta combine ####
if ( file.exists("04b_data_filt_and_combine/nmds_bray.RData") ) {
  load("04b_data_filt_and_combine/nmds_bray.RData")
  load("04b_data_filt_and_combine/nmds_wu.RData")
  
} else {
  set.seed(9234)
  nmds_bray <- isoMDS(as.dist(bray_filt), k=2)
  # nmds_uwu <- isoMDS(as.dist(uwu), k=2)
  nmds_wu <- isoMDS(as.dist(wu_filt), k=2)
  
  save(nmds_bray, file="04b_data_filt_and_combine/nmds_bray.RData")
  save(nmds_wu, file="04b_data_filt_and_combine/nmds_wu.RData")
  
}


meta_filt <- nmds_bray$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
  rename(bray_NMDS1=V1, bray_NMDS2=V2) %>%
  mutate(bray_stress = nmds_bray$stress) %>%
  full_join(meta_filt)
# 
# 
# meta_filt <- nmds_uwu$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   rename(uwu_NMDS1=V1, uwu_NMDS2=V2) %>%
#   mutate(uwu_stress = nmds_uwu$stress) %>%
#   full_join(meta_filt)

meta_filt <- nmds_wu$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
  rename(wu_NMDS1=V1, wu_NMDS2=V2) %>%
  mutate(wu_stress = nmds_wu$stress) %>%
  full_join(meta_filt)
gg_bray_nmds <- meta_filt %>% 
  ggplot() + geom_point(aes(x=bray_NMDS1, y=bray_NMDS2, col=HostPhylo), show.legend = FALSE)
gg_bray_nmds
ggsave("04b_data_filt_and_combine/NMDS_bray_all.png", 
       gg_bray_nmds
  )
# meta_filt %>% 
#   ggplot() + geom_point(aes(x=uwu_NMDS1, y=uwu_NMDS2, col=HostPhylo), show.legend = FALSE)
gg_wu_nmds <-  meta_filt %>% 
  ggplot() + geom_point(aes(x=wu_NMDS1, y=wu_NMDS2, col=HostPhylo), show.legend = FALSE)
gg_wu_nmds
ggsave("04b_data_filt_and_combine/NMDS_wu_all.png", 
      gg_wu_nmds
      )
## Weird streak
meta_filt %>% 
  filter(bray_NMDS1 > 0.75, bray_NMDS2 >0.5)
meta_filt %>% filter(HostPhylo=="Hyla_japonica") %>%
  # ggplot() + geom_point(aes(x=bray_NMDS1, y=bray_NMDS2, col=HostPhylo), show.legend = FALSE)
  ggplot() + geom_point(aes(x=wu_NMDS1, y=wu_NMDS2, col=HostPhylo), show.legend = FALSE)
# Remove this; looks super suspicious?


#### EDA pltos ####
## Look at richness across sample types
meta_filt %>%
  ggplot() + geom_violin(aes(x=Study, y=observed_features)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

meta_filt %>%
  ggplot() + geom_violin(aes(x=HostPhylo, y=observed_features)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

meta_filt %>%
  ggplot() + geom_violin(aes(x=HostPhylo, y=faith_pd)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Plot richness on tree

amphRichTips <- unlist(lapply(amphtree$tip.label, FUN=function(x) { mean((meta_filt[meta_filt$HostPhylo==x, "observed_features"]), na.rm=TRUE) }))
approxRich <- round((amphRichTips - min(amphRichTips, na.rm=TRUE))/max(amphRichTips, na.rm=TRUE) * 100)
redRamp <- colorRampPalette(c("grey","red"))
amphRichTipsCol <- redRamp(100)[approxRich]
png("04b_data_filt_and_combine/richness_amphibian_phylogeny.png", height=1000, width=1000)
plot(amphtree, tip.color = amphRichTipsCol)
dev.off()

meta_filt %>%
  ggplot() + geom_point(aes(x=as.numeric(Longitude), y=as.numeric(Latitude), col=log10(observed_features)), cex=4) 

meta_filt %>%
  ggplot() + geom_point(aes(x=as.numeric(Longitude), y=as.numeric(Latitude), col=HostPhylo), cex=4, show.legend = FALSE) 

####  Distribution of taxa ####
gg_distrSamples <- meta_filt %>%
  group_by(Host, Location) %>% summarize(nSamps=n()) %>% ungroup() %>%
  group_by(Host) %>% mutate(nSampHost = sum(nSamps)) %>% ungroup() %>%
  arrange(-nSampHost) %>% mutate(Host = factor(Host, levels=unique(Host))) %>%
  ggplot() + geom_bar(aes(x=Host, y=nSamps, fill=Location), stat="identity", show.legend = FALSE) +
  ylab("Number of Samples") + xlab("Host species") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_distrSamples
ggsave("04b_data_filt_and_combine/distribution_of_samples.png", height=5, width=10
       , gg_distrSamples)

# Filter otu and meta for prevalence testing
meta_cut <- meta_filt %>% filter(sampleid %in% colnames(otu))

# Combine OTU table in tsv form for FAPROTAX
otu_rare_withTaxa <- otu_rare
otu_rare_withTaxa$OTUID <-taxa_expand[match(rownames(otu_rare_withTaxa),taxa_expand$FeatureID),"Taxon"] %>% pull()
otu_rare_withTaxa <- otu_rare_withTaxa %>% select(OTUID, everything())

sink("04b_data_filt_and_combine/LOG.txt", append=TRUE)
print("Number of samples in meta_cut")
nrow(meta_cut)
print("Number of COMPLETE samples in meta_cut")
meta_cut %>% select(sampleid, Host, Latitude, Longitude, Month, Location) %>% drop_na() %>%
  nrow()
print("Number of Host")
meta_cut %>% select(Host) %>% distinct() %>% nrow()
print("Number of Studies")
meta_cut %>% select(Study) %>% distinct() %>% nrow()


sink()

###### Save #######
write.tree(tree, file="04b_data_filt_and_combine/downstream/tree_filt.nwk")
write.tree(tree_cal, file = "04b_data_filt_and_combine/downstream/tree_filt_cal.nwk")

write.table(otu_rare %>% rownames_to_column(var="ASVID"), file = "04b_data_filt_and_combine/downstream/otu_rare.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(otu %>% rownames_to_column(var="ASVID"), file = "04b_data_filt_and_combine/downstream/otu_cut.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(otu_rare_withTaxa, file = "04b_data_filt_and_combine/downstream/otu_rare_withTaxa.txt", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(meta_filt, file = "04b_data_filt_and_combine/downstream/meta_filt.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(meta_cut, file = "04b_data_filt_and_combine/downstream/meta_cut.txt", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(bray_filt%>%rownames_to_column(var="ROWNAMES"), file="04b_data_filt_and_combine/downstream/dm_bray.txt", row.names = FALSE, quote=FALSE, sep="\t")
write.table(wu_filt%>%rownames_to_column(var="ROWNAMES"), file="04b_data_filt_and_combine/downstream/dm_wu.txt", row.names = FALSE, quote=FALSE, sep="\t")
write.table(jac_filt%>%rownames_to_column(var="ROWNAMES"), file="04b_data_filt_and_combine/downstream/dm_jac.txt", row.names = FALSE, quote=FALSE, sep="\t")

write.table(taxa_expand, file="04b_data_filt_and_combine/downstream/taxonomy.txt", quote=FALSE, row.names=FALSE, sep="\t")
