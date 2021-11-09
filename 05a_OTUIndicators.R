#!bin/bash Rscript
library(MASS)
library(vegan)
library(ggfortify)
library(tidyverse)

dir.create("05a_OTUIndicators")
######## Loading #####
keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt") 
# filter(Host %in% keepAbund$Host& Location %in% keepAbund$Location, CollectionDate %in% keepAbund$CollectionDate)
otu <- read.delim("04b_data_filt_and_combine/downstream/otu_cut.txt", row.names=1) %>%
  select(one_of(meta$sampleid))
otu <- otu[rowSums(otu)>0,]

otuRA <- apply(otu, 2, function(x) x/sum(x))

#### otu indic ####
source("code/otu_indicator_functions.R")

if ( ! file.exists("05a_OTUIndicators/wScores1.RData") ) {
  otuRAFilt <- otuRA[apply(otuRA, 1, max)>0.01,]
  # otuRAFilt2 <- otuRA[apply(otuRA, 1, max)>0.05,]
  # TAKSE A LONG TIME
  wScores1 <- wCluster(otuRAFilt)
  # wScores0.5 <- wCluster(otuRAFilt2)
  save(wScores1, file="05a_OTUIndicators/wScores1.RData")
  # save(wScores0.5, file="05a_OTUIndicators/wScores0.5.RData")
} else {
  load("05a_OTUIndicators/wScores1.RData")
  # load("05a_OTUIndicators/wScores0.5.RData")
}

#### kmeans ####
#### Classify each OTU across sites as "absent", "low abundance" or "high abundance"
if ( ! file.exists("05a_OTUIndicators/otuClass_byotu.RData")) {
  otuRA_filt <- otuRA[rowSums(otuRA>0)>3,]
  otuClass_byotu <- matrix(nrow=dim(otuRA_filt)[1], ncol=dim(otuRA_filt)[2], dimnames = dimnames(otuRA_filt))
  for ( o in rownames(otuRA_filt) ) {
    # o =rownames(otuRA_filt)[1]
    kmeanstemp <- kmeans(log10(otuRA_filt[o,]+1), centers =3, nstart = 10)
    clusters=data.frame(clusters=kmeanstemp$cluster)
    clusterLegend =data.frame(clusters=order(kmeanstemp$centers), group=c(0,0.5,1)) 
    abund <- otuRA_filt[o,][rownames(clusters)]
    fullGroups <- cbind(clusters, abund) 
    fullGroups$group = clusterLegend[match(fullGroups$clusters, clusterLegend$clusters),"group"]
    # ggplot(fullGroups) + geom_histogram(aes(x=abund, fill=factor(group)))
    otuClass_byotu[o,] <- fullGroups[match(colnames(otuClass_byotu), rownames(fullGroups)), "group"]
  }
  save(otuClass_byotu, file="05a_OTUIndicators/otuClass_byotu.RData")
} else {
  load("05a_OTUIndicators/otuClass_byotu.RData")
}

# otuClass_nozeros <- otuClass_byotu[apply(otuClass_byotu, 1, max)>0,]
# otuNames <- rownames(otuClass_nozeros)
# comboMetaOTUClass <- t(otuClass_nozeros) %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   rename_at(vars(all_of(otuNames)), ~paste0("OTU_",otuNames)) %>%
#   right_join(meta) %>%
#   select(sampleid, Host, Location, biom_AnnualPrecip, biom_TempSeasonality, HabitatClass, Month, Latitude,Longitude, starts_with("OTU_")) %>%
#   pivot_longer(starts_with("OTU_"), names_to="OTU", values_to="Prev")
# comboMetaOTUClass %>%
#   filter(OTU%in% paste0("OTU_",otuNames[1:5])) %>%
#   group_by(OTU, Location) %>%
#   summarise(prev=sum(Prev>0)/n()) %>%
#   filter(prev>0) %>%
#   ggplot() + geom_point(aes(x=Location, y=OTU, cex=prev)) +
#   theme(axis.text.y = element_blank(), axis.text.x=element_text(angle=90))

## Do at sample level
if ( ! file.exists("05a_OTUIndicators/otuClass_bysamp.RData") ) {
  otuRA_filtSamp <- otuRA[,colSums(otuRA>0)>3]
  # hist(log10(colSums(otuRA_filtSamp>0)))
  otuClass_bysamp <- matrix(nrow=dim(otuRA_filtSamp)[1], ncol=dim(otuRA_filtSamp)[2], dimnames = dimnames(otuRA_filtSamp))
  for ( o in colnames(otuRA_filtSamp) ) {
    # o = colnames(otuRA_filt)[1]
    kmeanstemp <- kmeans(log10(otuRA_filtSamp[,o]+1), centers =3, nstart = 10)
    clusters=data.frame(clusters=kmeanstemp$cluster)
    clusterLegend =data.frame(clusters=order(kmeanstemp$centers), group=c(0,0.5,1)) 
    abund <- otuRA_filtSamp[,o][rownames(clusters)]
    fullGroups <- cbind(clusters, abund) 
    fullGroups$group = clusterLegend[match(fullGroups$clusters, clusterLegend$clusters),"group"]
    # View(fullGroups)
    # ggplot(fullGroups) + geom_histogram(aes(x=abund, fill=factor(clusters)))+xlim(0,0.05)
    otuClass_bysamp[,o] <- fullGroups[match(rownames(otuClass_bysamp), rownames(fullGroups)), "group"]
  }
  save(otuClass_bysamp, file="05a_OTUIndicators/otuClass_bysamp.RData")
} else {
  load("05a_OTUIndicators/otuClass_bysamp.RData")
}
# View(otuClass_bysamp)

# library(vegan)
# bc_all <- vegdist(t(otuRA_filt), method="bray")
# bc_sampclass <- vegdist(t(otuSampClass), method="bray")
# bc_class <- vegdist(t(otuClass), method="bray")
# 
# nmds_all <- isoMDS(bc_all, k=2)
# nmds_sampclass <- isoMDS(as.dist(bc_sampclass), k=2)
# nmds_class <- isoMDS(as.dist(bc_class), k=2)
# 
# nmds_all$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   left_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Location, pch=factor(Month)))+
#   facet_wrap(.~Host)
# 
# nmds_sampclass$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   left_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Host, pch=factor(Month)), cex=3) +
#   facet_wrap(.~Location)
# 
# nmds_all$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   left_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Host)) 
# nmds_class$points %>% as.data.frame() %>% rownames_to_column(var="sampleid") %>%
#   left_join(meta) %>%
#   ggplot() + geom_point(aes(x=V1, y=V2, col=Host)) 
# 
# tempAbundOTU <- names(sort(rowSums(otuRA), decreasing = TRUE))[1:10]
# hist(log10(otuRA[tempAbund[10],]))

