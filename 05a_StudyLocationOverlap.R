#!bin/bash

### Compare overlap in ASVs between studies and locations ###
library(ape)
library(lme4)
library(MASS, exclude = "select")
library(tidyverse)

dir.create("05a_StudyLocationOverlap")
#### Load
toKeep <- read.delim("04c_partition_data_for_downstream/multipleStudyLocList.txt")
meta <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt") %>% 
  inner_join(toKeep)
otu <- read.delim("04b_data_filt_and_combine/downstream/otu_cut.txt", row.names=1) %>%
  select(one_of(meta$sampleid))
otu <- otu[rowSums(otu)>0,]
# Get PA matrix
otuPA <- otu
otuPA[otuPA>0] <- 1
amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
bc <- read.delim("04b_data_filt_and_combine/downstream/dm_bray.txt", row.names = 1)[meta$sampleid, meta$sampleid]
wu <- read.delim("04b_data_filt_and_combine/downstream/dm_wu.txt", row.names = 1)[meta$sampleid, meta$sampleid]

### Look at samples ####
length(unique(meta$Host))
tree_filt <- keep.tip(amphtree, unique(meta$Host))
strHab <- meta$HabitatClass[match(tree_filt$tip.label, meta$Host)]
habcols <- data.frame(Habitat=c("Aquatic","Semi-aquatic","Arboreal","Terrestrial"), habCol=c("blue","green","seagreen","darkred"))
strHabCol <- left_join(data.frame(Habitat=strHab),habcols) %>% pull(habCol)
png("05a_StudyLocationOverlap/amphibians_included_tree.png", height=1000, width = 800)
plot(tree_filt, tip.color=strHabCol)
dev.off()
# 
# nmds_wu <- isoMDS(as.dist(wu), k=2)
# nmds_bc <- isoMDS(as.dist(bc), k=2)
# 
# nmds_wu$points %>% as.data.frame() %>%
#   rownames_to_column(var="sampleid") %>% full_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Location), cex=3) 
# nmds_wu$points %>% as.data.frame() %>%
#   rownames_to_column(var="sampleid") %>% full_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Host), cex=3) 
# 
# 
# nmds_bc$points %>% as.data.frame() %>%
#   rownames_to_column(var="sampleid") %>% full_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Location), cex=3) 
# nmds_bc$points %>% as.data.frame() %>%
#   rownames_to_column(var="sampleid") %>% full_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Host), cex=3) 
# 

### Shared ASVs between ALL samples
# nSharedMat_all <- t(as.matrix(otuPA))%*%as.matrix(otuPA)
# nTotalMat_all <- nrow(as.matrix(otuPA)) - t(1-as.matrix(otuPA))%*%(1-as.matrix(otuPA))

### Get shared ASVs for every SpxLocxDate combo; then can aggregate within- and between- Study after

####  Shared ASVs between samples within SpxLocxDate  ####
meta_adj <- meta %>% unite(Study,Location, CollectionDate, col="StLocDate", remove=FALSE) 

allSp <- unique(meta_adj$Host)
allMatCount <- list()
for ( sp in allSp ) {
  meta_temp <- meta_adj %>% filter(Host==sp) %>% select(sampleid, Study, Location, Host, StLocDate)
  allMatCount[[sp]] <- list()
  # Combined otu with meta?
  otuPA_temp <- otuPA[,meta_temp$sampleid]
  otuPA_mat <- otuPA_temp[rowSums(otuPA_temp)>0,]
  #Number of shared ASV
  nSharedmat <- t(as.matrix(otuPA_mat))%*%as.matrix(otuPA_mat)
  # Turn 0to1 and vv; then do same thing to calculate total shared ASVs
  nTotalMat <- nrow(otuPA_mat)- t(1-as.matrix(otuPA_mat))%*%(1-as.matrix(otuPA_mat))
  # Check this
  # nTotalMat["SRR3502315","SRR7614999"]
  # sum(rowSums(otuPA_mat[,c("SRR3502315","SRR7614999")])>0)
  
  allMatCount[[sp]][["SharedASVs"]] <- nSharedmat
  allMatCount[[sp]][["TotalASVs"]] <- nTotalMat
  
}



# Combine into long plots for all species
sharedASVs_long <- data.frame()
for (sp in allSp) {
  tempMat <- allMatCount[[sp]]
  nASVs <- tempMat$SharedASVs
  tASVs <- tempMat$TotalASVs
  # pASVs <- nASVs/tASVs
  nASVs[lower.tri(tempMat$SharedASVs)] <- NA
  pairs <- t(combn(rownames(nASVs), 2))
  tempLong <- data.frame(Host=sp, pairs, nASV=nASVs[pairs], tASV=tASVs[pairs]) 
  if (nrow(sharedASVs_long)==0) {
    sharedASVs_long <- tempLong
  } else {
    sharedASVs_long <- rbind(sharedASVs_long, tempLong)
  }
}

sharedASVs_long_adj <- sharedASVs_long %>%
  mutate(Loc1=meta_adj[match(X1, meta_adj$sampleid),"Location"]
         , Loc2=meta_adj[match(X2, meta_adj$sampleid),"Location"]
         , Date1=meta_adj[match(X1, meta_adj$sampleid),"CollectionDate"]
         , Date2=meta_adj[match(X2, meta_adj$sampleid),"CollectionDate"]
         , Study1=meta_adj[match(X1, meta_adj$sampleid),"Study"]
         , Study2=meta_adj[match(X2, meta_adj$sampleid),"Study"]
  ) %>%
  mutate(Date1=ifelse(is.na(Date1), "2021", Date1), Date2=ifelse(is.na(Date2), "2021", Date2)) %>%
  mutate(Location=Loc1==Loc2, Date=Date1==Date2, Study=Study1==Study2) %>%
  mutate(pASV=nASV/tASV, logpASV=log(nASV+1)/log(tASV+1)) %>%
  mutate(SampleBout = Location*Date)


sharedASVs_long_adj %>%
  ggplot() + geom_violin(aes(x=Host, y=logpASV)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

sharedASVs_long_adj %>%
  group_by(Host) %>%
  mutate(meantASV=mean(tASV, na.rm=TRUE)) %>% ungroup() %>%
  ggplot() + geom_violin(aes(x=Host, y=log10(nASV+1))) +
  geom_point(aes(x=Host, y=log10(meantASV+1)), col="red")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


sharedASVs_long_adj %>%
  group_by(Host) %>%
  mutate(meanlogtASV=mean(logpASV, na.rm=TRUE)) %>% ungroup() %>%
  arrange(meanlogtASV) %>%
  mutate(Host=factor(Host, levels=unique(Host))) %>%
  ggplot() + geom_boxplot(aes(x=Host, y=logpASV)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

### Correlating with enviro var
allbiom <- meta_adj %>% select(starts_with("biom_")) %>% colnames()
sharedASVs_long_adj %>%
  filter(Location, Date) %>% # Only comparing within-location
  mutate(Location=Loc1, Date=Date1) %>% select(Host, Location, Date,nASV, tASV, pASV, logpASV) %>%
  left_join(meta_adj)  %>%
  filter(!is.na(biom_AnnualMeanTemp)) %>%
  group_by(Host,Location,Date) %>%
  summarize(meanlogpASV=mean(logpASV), sdlogpASV=sd(logpASV) ) %>%
  left_join(meta%>%select(Location, starts_with("biom_"))) %>%
  ggplot() + geom_point(aes(x=biom_TempSeasonality, y=meanlogpASV, col=Host), show.legend = FALSE) +
  geom_smooth(aes(x=biom_TempSeasonality, y=meanlogpASV), method="lm") +
  geom_errorbar(aes(x=biom_TempSeasonality, ymin=meanlogpASV-sdlogpASV, ymax=meanlogpASV+sdlogpASV))
allbiom

meta%>%select(Location, starts_with("biom_")) %>% tibble()
# View(sharedASVs_long)
sharedASVs_long_adj %>%
  ggplot() + geom_violin(aes(x=Study, y=pASV)) +
  facet_wrap(.~Host)

lm_samplebyhost <- lm(logpASV ~ Host*SampleBout, data=sharedASVs_long_adj)
sum_samplebyhost <- summary(lm_samplebyhost)$coefficients %>% as.data.frame() %>%
  rownames_to_column(var="Predictor") %>%
  rename_at(vars("Estimate", "Std. Error", "t value", "Pr(>|t|)"), ~c("b","se","t","p")) %>%
  separate(Predictor, into=c("Host","SampleBout"), sep=":", remove=FALSE) %>%
  mutate(SampleBout = ifelse(Host=="SampleBout", "SampleBout", SampleBout), Host=gsub("Host","",Host)) %>%
  mutate(Host= ifelse(Host%in% c("SampleBout", "(Intercept)"), "Intercept",Host) 
         , SampleBout = ifelse(is.na(SampleBout), "EffectOfHost", "EffectOfSampleBout")) %>%
  rename(Type=SampleBout)
#
sum_samplebyhost %>% View()

#Look at tree
greenRamp <- colorRampPalette(c("red","green"))
greenRampCol <- greenRamp(100)
blueRamp <- colorRampPalette(c("yellow","blue"))
blueRampCol <- blueRamp(100)
## EFFECT OF HOST ON SHARED ASVs
b_hostEff <- sum_samplebyhost[match(tree_filt$tip.label, sum_samplebyhost %>% filter(Type=="EffectOfHost") %>% pull(Host)), "b"]
b_hostEff_st <- round((b_hostEff-min(b_hostEff, na.rm = TRUE))/(max(b_hostEff,na.rm=TRUE) -min(b_hostEff, na.rm = TRUE))*100)
b_hostEff_col <- greenRampCol[b_hostEff_st]
plot(tree_filt, tip.color = b_hostEff_col)

b_sampleEff <- sum_samplebyhost[match(tree_filt$tip.label, sum_samplebyhost %>% filter(Type=="EffectOfSampleBout") %>% pull(Host)), "b"]
b_sampleEff_st <- round((b_sampleEff-min(b_sampleEff, na.rm = TRUE))/(max(b_sampleEff,na.rm=TRUE) -min(b_sampleEff, na.rm = TRUE))*100)
b_sampleEff_col <- blueRampCol[b_sampleEff_st]
plot(tree_filt, tip.color = b_sampleEff_col)



##### Prevalence vs correlation ####
otu_RA <- apply(otu, 2, function(x) x/sum(x))
meta_adj %>% mutate(Intercept=1) %>% select(Intercept)

otu_RA[1,]

