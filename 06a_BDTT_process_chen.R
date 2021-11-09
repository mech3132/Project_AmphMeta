#!bin/bash

#### Control on BDTT-- used back cut tips


library(ape)
library(castor)
library(caper)
library(tidyverse)
dir.create("06a_BDTT_process_chen")

### Load ###
# keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt")
taxa <- read.delim("04b_data_filt_and_combine/downstream/taxonomy.txt")
bacttree <- read.tree("04b_data_filt_and_combine/downstream/tree_filt.nwk")
bacttree_bdtt <- read.tree("04b_data_filt_and_combine/downstream/tree_filt_cal.nwk")

load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_pooled_chen/allDat.RData")
allDat_chen <- allDat
gg_full_chen <- allDat_chen %>%
  filter(type=="Full") %>%
  ggplot(aes(x=slice, y=abs(Coef), col=predictor)) +
  geom_point() + geom_line() +
  # facet_wrap(.~type)+ 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  xlim(0,1)
gg_full_chen
ggsave(filename = "06a_BDTT_process_chen/correlations_pooled_chen_coef.png"
       ,gg_full_chen)



# 
# load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_rare_chen/allDat.RData")
# allDat_chen <- allDat
# # load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen/allDat.RData")
# # allDat_bdtt <- allDat
# 
# 
# 
# ##### Chen #####
# 
# gg_cordat_bytype_chen <- allDat_chen %>%
#   ggplot(aes(x=slice, y=R2, col=group)) +
#   geom_point() + geom_line() +
#   facet_wrap(.~type)+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
#   xlim(0,1)
# gg_cordat_bytype_chen
# ggsave(filename = "06a_BDTT_process_chen/correlations_bytype_bray_rare_chen.png"
#        ,gg_cordat_bytype_chen)
# 
# gg_cordat_bygroup_chen <- allDat_chen %>%
#   ggplot(aes(x=slice, y=R2, col=type)) +
#   geom_point() + geom_line() +
#   facet_wrap(.~group)+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
#   xlim(0,1)
# gg_cordat_bygroup_chen
# ggsave(filename = "06a_BDTT_process_chen/correlations_bygroup_bray_rare_chen.png"
#        ,gg_cordat_bygroup_chen)
# 
# 
# # Find maxes
# maxR2_chen <- allDat_chen %>% group_by(type, group) %>%
#   mutate(maxR2=max(R2), thisMax=maxR2==R2) %>%
#   filter(thisMax)
# 
# maxR2_chen_filt <- maxR2_chen %>% filter(group=="full") %>% distinct() %>% arrange(-R2)
# write.table(maxR2_chen_filt, file="06a_BDTT_process_chen/table_maxR2_bdtt.txt", quote=FALSE, row.names = FALSE, sep="\t")
# 
# 
# ### Find out depth of geo collapse)
# slices=round(seq(from=0, to=1.5, length.out=40),2)
# if ( !file.exists("06a_BDTT_process_chen/allSlicesStats.RData") ) {
#   allSlicesStats <- data.frame(matrix(ncol=4, nrow=0, dimnames = list(NULL, c("percid","taxlevel","ncollapsed","slice"))))
#   for (r in slices[1:10]) {
#     print(r)
#     collapsedTree <- collapse_tree_at_resolution(bacttree, resolution = r, rename_collapsed_nodes = TRUE)
#     nodesToCollapse <- collapsedTree$collapsed_nodes+length(bacttree$tip.label)
#     names(nodesToCollapse) <- bacttree$tip.label[collapsedTree$farthest_tips]
#     simBetweenTips <- c()
#     simTaxLevel <- c()
#     ncollapsed <- c()
#     # length(nodesToCollapse)
#     for ( n in names(nodesToCollapse)) {
#       # n <- names(nodesToCollapse)[1]
#       allMembers <- clade.members(nodesToCollapse[[n]], phy=bacttree, tip.labels = TRUE)
#       ncollapsed <- c(ncollapsed, length(allMembers))
#       # Get common taxonomic level
#       temptax <- t(taxa[taxa$FeatureID%in%allMembers,c("Domain","Phylum","Class","Order","Family","Genus","Species")])
#       simTaxLevel <- c(simTaxLevel, max(which(apply(temptax, 1, function(x) length(unique(x)))==1)))
#       # Get percent identity
#       strSeq <- strsplit(allMembers, split="")
#       df <- data.frame(matrix(unlist(strSeq), ncol=length(strSeq), byrow=FALSE))
#       percid <- sum(apply(df, 1, function(x) length(unique(x)))==1)/nrow(df)
#       simBetweenTips <- c(simBetweenTips, percid)
#     }
#     allTemp <- data.frame(percid=simBetweenTips, taxlevel=simTaxLevel, ncollapsed=ncollapsed, slice=r)
#     allSlicesStats <- rbind(allSlicesStats, allTemp)
#   }
#   save(allSlicesStats, file="06a_BDTT_process_chen/allSlicesStats.RData")
# } else {
#   load("06a_BDTT_process_chen/allSlicesStats.RData")
# }
# 
# 
# ggsave(filename = "06a_BDTT_process_chen/percent_dff_btwn_tips_0.27slice.png",
#        allSlicesStats %>%
#          filter(slice==0.27) %>%
#          ggplot() + geom_histogram(aes(x=percid))+
#          xlab("Percent identity between collapsed tips")
# )
# 
# sumTab_chen <- allSlicesStats %>%
#   group_by(slice) %>%
#   summarise(meanpercid=mean(percid), medpercid=median(percid), medtaxlevel=median(taxlevel), meanncollapsed=mean(ncollapsed))
# write.table(sumTab_chen, file = "06a_BDTT_process_chen/table_collapsedeptsh_and_values_chen.txt",quote = FALSE, row.names = FALSE, sep="\t")
# 
# 
# ggslicepercid <- allSlicesStats %>%
#   # filter(percid>0.3) %>%
#   ggplot() + geom_boxplot(aes(x=factor(slice), y=percid)) +
#   xlab("Slice depth into phylogenetic tree") +
#   ylab("Percent identity between collapsed tips")
# ggslicepercid
# ggsave(filename = "06a_BDTT_process_chen/percent_identity_across_slices.png"
#        ,ggslicepercid)
# 
# ggntipscoll <-  allSlicesStats %>%
#   ggplot() + geom_boxplot(aes(x=factor(slice), y=ncollapsed)) +
#   xlab("Slice depth into phylogenetic tree") +
#   ylab("Number taxa collapsed per node\n(log10 scale)") +
#   scale_y_log10()
# ggntipscoll
# ggsave(filename = "06a_BDTT_process_chen/number_tips_collapsed.png"
#        ,ggntipscoll)
# 
# 
# gg_taxcollapse <- allSlicesStats %>%
#   # filter(percid>0.3) %>%
#   mutate(taxlevel=ifelse(taxlevel==1,"Kingdom",ifelse(taxlevel==2, "Phylum"
#                                                       , ifelse(taxlevel==3,"Class"
#                                                                , ifelse(taxlevel==4,"Order"
#                                                                         , ifelse(taxlevel==5,"Family"
#                                                                                  ,ifelse(taxlevel==6,"Genus","Species"))))))) %>%
#   mutate(taxlevel=factor(taxlevel, levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
#   ggplot() + geom_bar(aes(x=factor(slice), fill=factor(taxlevel))) +
#   xlab("Slice depth into phylogenetic tree") +
#   ylab("Number of nodes")+
#   labs(fill="Highest resolution in taxonomy\nshared by all collapsed nodes")
# ggsave(filename = "06a_BDTT_process_chen/taxonomy_collapsed_byslice.png"
#        ,gg_taxcollapse)
# 
# # MODEsimilarity <- names(which(table(simBetweenTips)==max(table(simBetweenTips))))
# # MEDsimilarity <- median(simBetweenTips)
# # MODEsimilarity
# # MEDsimilarity
# 
# #### Shuffle ####
# load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_rare_chen_SHUFFLE/CorsbcphyloDist_Full.RData")
# Cors_phylo_full <- Cors_full
# load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_rare_chen_SHUFFLE/CorsbcecoDist_Full.RData")
# Cors_eco_full <- Cors_full
# load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_rare_chen_SHUFFLE/CorsbcgeoDist_Full.RData")
# Cors_geo_full <- Cors_full
# 
# colnames(Cors_phylo_full) <- slices
# phylocorr <- Cors_phylo_full %>% t() %>% as.data.frame() %>% rownames_to_column("slice") %>%
#   mutate(slice=as.numeric(slice), type="phyloDist", group="shuffleFull")
# colnames(Cors_eco_full) <- slices
# ecocorr <- Cors_eco_full %>% t() %>% as.data.frame() %>% rownames_to_column("slice") %>%
#   mutate(slice=as.numeric(slice), type="ecoDist", group="shuffleFull")
# colnames(Cors_geo_full) <- slices
# geocorr <- Cors_geo_full %>% t() %>% as.data.frame() %>% rownames_to_column("slice") %>%
#   mutate(slice=as.numeric(slice), type="geoDist", group="shuffleFull")
# 
# ggControlShuffle <- rbind(phylocorr, ecocorr, geocorr) %>% full_join(allDat) %>%
#   filter(group %in% c("full","shuffleFull"), type %in% c("phyloDist","ecoDist","geoDist")) %>%
#   mutate(type=ifelse(type=="ecoDist","Correlation with\nHabitat Ecology",
#                      ifelse(type=="geoDist","Correlation with\ngeographic distance",
#                             ifelse(type=="phyloDist","Correlation with\nhost phylogeny",NA)))) %>%
#   mutate(group=ifelse(group=="full","Full structured data", "Site IDs shuffled")) %>%
# ggplot() + geom_line(aes(x=slice, y=R2)) +
#   facet_grid(group~type, scales="free") +
#   xlab("Depth of slice into phylogenetic tree")
# ggControlShuffle
# ggsave(filename = "06a_BDTT_process_chen/shuffled_controls.png"
#        ,ggControlShuffle)

