#!bin/bash
library(ape)
library(castor)
library(caper)
library(tidyverse)
dir.create("06a_BDTT_process")

### Load ###
# keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt")
taxa <- read.delim("04b_data_filt_and_combine/downstream/taxonomy.txt")
bacttree <- read.tree("04b_data_filt_and_combine/downstream/tree_filt.nwk")
bacttree_bdtt <- read.tree("04b_data_filt_and_combine/downstream/tree_filt_cal.nwk")

# load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen/allDat.RData")
# allDat_bdtt <- allDat
# Load each individually
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_climFull/allDat.RData")
climDat <- allDat
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_phyloFull/allDat.RData")
phyloDat <- allDat
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_ecoFull/allDat.RData")
ecoDat <- allDat
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_geoFull/allDat.RData")
geoDat <- allDat
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_timeFull/allDat.RData")
timeDat <- allDat
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_chen_studyFull/allDat.RData")
studyDat <- allDat

allDat_bdtt <- rbind(phyloDat, ecoDat, geoDat, timeDat, climDat, studyDat) %>% 
  filter(predictor!="Int")
##### BDTT version ######

gg_cordat_separateR2_bdtt <- allDat_bdtt %>% 
  ggplot(aes(x=slice, y=R2, col=predictor)) +
  geom_point(aes(pch=pval_coef<0.05)) + geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  xlim(0,1)+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_cordat_separateR2_bdtt
ggsave(filename = "06a_BDTT_process/correlations_separate_bdtt_R2.png"
       ,gg_cordat_separateR2_bdtt)

gg_cordat_separatecoef_bdtt <- allDat_bdtt %>% 
  ggplot(aes(x=slice, y=Coef, col=predictor)) +
  geom_point(aes(pch=pval_coef<0.05)) + geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  xlim(0,1)+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_cordat_separatecoef_bdtt
ggsave(filename = "06a_BDTT_process/correlations_separate_bdtt_coef.png"
       ,gg_cordat_separatecoef_bdtt)

# Find maxes
maxR2_bdtt <- allDat_bdtt %>% group_by(type, predictor) %>%
  mutate(maxR2=max(R2), thisMax=maxR2==R2) %>%
  filter(thisMax) 

maxR2_bdtt_filt <- maxR2_bdtt %>% distinct() %>% arrange(-R2)
write.table(maxR2_bdtt_filt, file="06a_BDTT_process/table_maxR2_bdtt.txt", quote=FALSE, row.names = FALSE, sep="\t")
maxR2_bdtt_filt

save(allDat_bdtt, file = "06a_BDTT_process/allDat_bdtt_separate.RData")
### BDTT version
slices_bdtt <- round(seq(0,2,length.out=40),2)
source("code/BDTT-master/BDTT_functions_MYC.R")
if ( !file.exists("06a_BDTT_process/allSlicesStats_bdtt.RData") ) {
  allSlicesStats_bdtt <- data.frame(matrix(ncol=4, nrow=0, dimnames = list(NULL, c("percid","taxlevel","ncollapsed","slice"))))
  for (r in slices_bdtt[1:10]) {
    print(r)
    edgeInfo=getTreeInfo(bacttree_bdtt)
    brCorres=GetBranchNode(slice=r,edgeInfo=edgeInfo)
    nodesToCollapse=brCorres[,2]    
    #print("getting OTUS host matrix")
    # OTUnames=as.character(brCorres[,2])
    
    simBetweenTips <- c()
    simTaxLevel <- c()
    ncollapsed <- c()
    # length(nodesToCollapse)
    for ( n in 1:length(nodesToCollapse)) {
      # n <- 2
      allMembers <- clade.members(nodesToCollapse[n], phy=bacttree_bdtt, tip.labels = TRUE)
      ncollapsed <- c(ncollapsed, length(allMembers))
      # Get common taxonomic level
      temptax <- t(taxa[taxa$FeatureID%in%allMembers,c("Domain","Phylum","Class","Order","Family","Genus","Species")])
      simTaxLevel <- c(simTaxLevel, max(which(apply(temptax, 1, function(x) length(unique(x)))==1)))
      # Get percent identity
      strSeq <- strsplit(toupper(allMembers), split="")
      df <- data.frame(matrix(unlist(strSeq), ncol=length(strSeq), byrow=FALSE))
      percid <- sum(apply(df, 1, function(x) length(unique(x)))==1)/nrow(df)
      simBetweenTips <- c(simBetweenTips, percid)
    }
    allTemp <- data.frame(percid=simBetweenTips, taxlevel=simTaxLevel, ncollapsed=ncollapsed, slice=r)
    allSlicesStats_bdtt <- rbind(allSlicesStats_bdtt, allTemp)
  }
  save(allSlicesStats_bdtt, file="06a_BDTT_process/allSlicesStats_bdtt.RData")
} else {
  load("06a_BDTT_process/allSlicesStats_bdtt.RData")
}

ggsave(filename = "06a_BDTT_process/percent_dff_btwn_tips_0.27slice_bdtt.png", 
       allSlicesStats_bdtt %>%
         filter(slice==slices_bdtt[7]) %>%
         ggplot() + geom_histogram(aes(x=percid))+
         xlab("Percent identity between collapsed tips")
)

sumTab <- allSlicesStats_bdtt %>%
  group_by(slice) %>%
  summarise(meanpercid=mean(percid), medpercid=median(percid), medtaxlevel=median(taxlevel), meanncollapsed=mean(ncollapsed))
write.table(sumTab, file = "06a_BDTT_process/table_collapsedeptsh_and_values_bdtt.txt",quote = FALSE, row.names = FALSE, sep="\t")

ggslicepercid_bdtt <- allSlicesStats_bdtt %>%
  # filter(percid>0.3) %>%
  ggplot() + geom_boxplot(aes(x=factor(slice), y=percid)) +
  xlab("Slice depth into phylogenetic tree") +
  ylab("Percent identity between collapsed tips")
ggslicepercid_bdtt
ggsave(filename = "06a_BDTT_process/percent_identity_across_slices_bdtt.png"
       ,ggslicepercid_bdtt)

ggntipscoll_bdtt <-  allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>%
  ggplot() + geom_boxplot(aes(x=factor(slice), y=ncollapsed)) +
  xlab("Slice depth into phylogenetic tree") +
  ylab("Number taxa collapsed per node\n(log10 scale)") + 
  scale_y_log10()
ggntipscoll_bdtt
ggsave(filename = "06a_BDTT_process/number_tips_collapsed_bdtt.png"
       ,ggntipscoll_bdtt)


gg_taxcollapse_bdtt <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>%
  mutate(taxlevel=ifelse(taxlevel==1,"Kingdom",ifelse(taxlevel==2, "Phylum"
                                                      , ifelse(taxlevel==3,"Class"
                                                               , ifelse(taxlevel==4,"Order"
                                                                        , ifelse(taxlevel==5,"Family"
                                                                                 ,ifelse(taxlevel==6,"Genus","Species"))))))) %>%
  mutate(taxlevel=factor(taxlevel, levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot() + geom_bar(aes(x=factor(slice), fill=factor(taxlevel))) +
  xlab("Slice depth into phylogenetic tree") +
  ylab("Number of nodes")+
  labs(fill="Highest resolution in taxonomy\nshared by all collapsed nodes")
gg_taxcollapse_bdtt
ggsave(filename = "06a_BDTT_process/taxonomy_collapsed_byslice_bdtt.png"
       ,gg_taxcollapse_bdtt)

# MODEsimilarity <- names(which(table(simBetweenTips)==max(table(simBetweenTips))))
# MEDsimilarity <- median(simBetweenTips)
# MODEsimilarity
# MEDsimilarity


### Look at pooled ####
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_ALLTOGETHER/allDat.RData")
allDat_pooled_bdtt <- allDat %>% mutate(drop="fulldataset")
gg_allpred_together <- allDat_pooled_bdtt %>% filter(predictor!="Int") %>%
  ggplot() + geom_line(aes(x=slice, y=abs(Coef), col=predictor))+
  geom_point(aes(x=slice, y=abs(Coef), col=predictor, pch=pval_coef<0.05))+
  facet_wrap(.~type) +
  xlim(0,1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_allpred_together
ggsave(filename = "06a_BDTT_process/correlations_allpooledpredictors.png"
       ,gg_allpred_together)

allDat_pooled_bdtt %>% filter(predictor!="Int") %>%
  select(slice, R2, type,pval) %>% distinct() %>%
  ggplot() + geom_line(aes(x=slice, y=R2, col=type))+
  geom_point(aes(x=slice, y=R2, col=type, pch=pval<0.05))+
  xlim(0,1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))

### Drop systematically ####
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_NoPhylo/allDat.RData")
allDat_nophylo <- allDat %>% mutate(drop="nophylo")
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_NoGeo/allDat.RData")
allDat_nogeo <- allDat%>% mutate(drop="nogeo")
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_NoClim/allDat.RData")
allDat_noclim <- allDat%>% mutate(drop="noclim")
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_NoEco/allDat.RData")
allDat_noeco <- allDat%>% mutate(drop="noeco")
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_Notime/allDat.RData")
allDat_notime <- allDat%>% mutate(drop="notime")
load("05b_BDTT_slicing_originalBDTT/Bray_rare_chen/Correlations_bray_rare_NoStudy/allDat.RData")
allDat_nostudy <- allDat%>% mutate(drop="nostudy")

allDat_pooled_withdrop <- rbind(allDat_pooled_bdtt,allDat_nophylo, allDat_nogeo, allDat_noclim, allDat_noeco, allDat_notime, allDat_nostudy) %>%
  filter(predictor!="Int", type=="Full")
View(allDat_pooled_withdrop)


gg_allpooled <- allDat_pooled_withdrop %>%
  filter(drop=="fulldataset") %>%
  ggplot() + geom_line(aes(x=slice, y=abs(Coef), col=predictor))+
  geom_point(aes(x=slice, y=abs(Coef), col=predictor, pch=pval_coef<0.05))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  xlim(0,1) +
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_allpooled
ggsave(filename = "06a_BDTT_process/correlations_pooledpredictors_coef.png", height=5, width=8
       ,gg_allpooled)



gg_droptest <- allDat_pooled_withdrop %>%
  ggplot() + geom_line(aes(x=slice, y=abs(Coef), col=drop))+
  geom_point(aes(x=slice, y=abs(Coef), col=drop, pch=pval_coef<0.05))+
  facet_wrap(.~predictor) +
  xlim(0,1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_droptest
ggsave(filename = "06a_BDTT_process/correlations_bypredictor_droptest.png", height=5, width=8
       ,gg_droptest)


gg_droptest_bydrop <- allDat_pooled_withdrop %>%
  ggplot() + geom_line(aes(x=slice, y=abs(Coef), col=predictor))+
  geom_point(aes(x=slice, y=abs(Coef), col=predictor, pch=pval_coef<0.05))+
  facet_wrap(.~drop) +
  xlim(0,1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_droptest_bydrop
ggsave(filename = "06a_BDTT_process/correlations_bydrop_droptest.png", height=5, width=8
       ,gg_droptest_bydrop)


gg_droptest_R2 <- allDat_pooled_withdrop %>% 
  select(slice,R2,drop, pval) %>%
  ggplot() + geom_line(aes(x=slice, y=R2, col=drop))+
  geom_point(aes(x=slice, y=R2, col=drop, pch=pval<0.05))+
  xlim(0,1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_shape_manual(values=c("TRUE"=19,"FALSE"=21))
gg_droptest_R2
ggsave(filename = "06a_BDTT_process/correlations_R2r_droptest.png"
       ,gg_droptest_R2)



full_tab <- allDat_pooled_bdtt %>% filter(predictor!="Int", type=="Full") %>%
  select(slice, R2) %>% distinct() 
full_tab
write.table(full_tab, file="06a_BDTT_process/tab_combinedMRM_R2.txt", quote=FALSE, row.names=FALSE, sep="\t")

full_tab_with_drops <- allDat_pooled_withdrop %>% 
  select(slice,R2,drop, pval) %>% distinct() %>%
  filter(slice==0.05)%>% arrange(R2)
full_tab_with_drops
write.table(full_tab_with_drops, file="06a_BDTT_process/tab_combinedMRM_withdrops_R2.txt", quote=FALSE, row.names=FALSE, sep="\t")

full_tab_with_drops_25 <- allDat_pooled_withdrop %>% 
  select(slice,R2,drop, pval) %>% distinct() %>%
  filter(slice==0.26)%>% arrange(R2)
full_tab_with_drops_25
write.table(full_tab_with_drops_25, file="06a_BDTT_process/tab_combinedMRM_withdrops_at0.26_R2.txt", quote=FALSE, row.names=FALSE, sep="\t")

## Get R2 explained by each component
R2Full <- allDat_pooled_withdrop %>% filter(drop=="fulldataset") %>%
  select(slice,R2) %>% distinct() %>% rename(fullR2=R2)

changeinR2_over_slice_perpred <- allDat_pooled_withdrop  %>% filter(drop!="fulldataset") %>%
  left_join(R2Full) %>%
  select(R2, slice, drop, fullR2) %>% distinct() %>%
  mutate(changeinR2=fullR2-R2) 
changeinR2_over_slice_perpred
write.table(changeinR2_over_slice_perpred, file="06a_BDTT_process/tab_changeinR2_droptest.txt", quote=FALSE, row.names=FALSE, sep="\t")

gg_changeinR2 <- changeinR2_over_slice_perpred %>%
  rename(Predictor=drop) %>%
  mutate(Predictor=gsub("no","",Predictor)) %>%
  ggplot() + geom_line(aes(x=slice, y=changeinR2, col=Predictor)) +
  facet_wrap(.~Predictor, scales = "free")

gg_changeinR2 
ggsave(filename = "06a_BDTT_process/change_in_R2_over_slice.png", height=4, width=8
       ,gg_changeinR2)
# scaled
gg_changeinR2_scaled <- changeinR2_over_slice_perpred %>%
  rename(Predictor=drop) %>%
  mutate(Predictor=gsub("no","",Predictor)) %>%
  ggplot() + geom_line(aes(x=slice, y=changeinR2, col=Predictor)) +
  # facet_wrap(.~Predictor) +
  scale_y_log10() +ylab("Change in R2, log10 scale")

gg_changeinR2_scaled 
ggsave(filename = "06a_BDTT_process/change_in_R2_over_slice_scaled.png", height=4, width=8
       ,gg_changeinR2_scaled)


save(allDat_pooled_withdrop, file = "06a_BDTT_process/allDat_pooled_withdrop.RData")