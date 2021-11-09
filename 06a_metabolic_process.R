#!bin/bash
library(ecodist)
library(tidyverse)
##### Correlating functions now ######
dir.create("06a_metabolic_process")

# load("05a_metagenomic_inference/allMetaCor.RData")
load("05a_metagenomic_inference/pca_fapro_filt_relabund.RData")
load("05a_metagenomic_inference/pca_ec_filt_relabund.RData")
load("05a_metagenomic_inference/pca_ko_filt_relabund.RData")
load("05a_metagenomic_inference/allMetaCor_bc.RData")
load("05a_metagenomic_inference/allMetaCor_pooled_bc.RData")
load("05a_metagenomic_inference/allMetaCor_dropone.RData")

load("04d_host_and_enviro_matrices/RData/phyloDist.RData")
phyloDist <- as.matrix(phyloDist)
load("04d_host_and_enviro_matrices/RData/geoDist.RData")
geoDist <- as.matrix(geoDist)
load("04d_host_and_enviro_matrices/RData/ecoDist.RData")
ecoDist <- as.matrix(ecoDist)
load("04d_host_and_enviro_matrices/RData/climDist.RData")
climDist <- as.matrix(climDist)
load("04d_host_and_enviro_matrices/RData/timeDist.RData")
timeDist <- as.matrix(timeDist)



gg_singlecorrelations_R2 <- allMetaCor_bc %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef, n=5)) %>%
  filter(metabType!="fapro") %>%
  ggplot() + geom_bar(aes(x=type, y=R2, fill=type, col=padj<0.05), stat="identity" )+
  facet_grid(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_color_manual(values=c("FALSE"="white", "TRUE"="black"))
gg_singlecorrelations_R2
ggsave("06a_metabolic_process/plot_singleMRMcorrelations_R2.png"
       ,gg_singlecorrelations_R2)


gg_singlecorrelations_coef <- allMetaCor_bc %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef, n=5)) %>%
  filter(metabType!="fapro") %>%
  ggplot() + geom_bar(aes(x=type, y=coef, fill=type, col=padj<0.05), stat="identity" )+
  facet_grid(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_color_manual(values=c("FALSE"="white", "TRUE"="black"))
gg_singlecorrelations_coef
ggsave("06a_metabolic_process/plot_singleMRMcorrelations_coef.png"
       ,gg_singlecorrelations_coef)

gg_pooledcorrelations_coef <- allMetaCor_pooled_bc %>%
  filter(metabType!="fapro") %>%
  ggplot() + geom_bar(aes(x=predictor, y=abs(coef), fill=predictor, col=pval_coef<0.05), stat="identity" )+
  facet_grid(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_color_manual(values=c("FALSE"="white", "TRUE"="black"))
gg_pooledcorrelations_coef
ggsave("06a_metabolic_process/plot_pooledMRMcorrelations_coef.png"
       ,gg_pooledcorrelations_coef)

### DROP ONE:
gg_R2drop_meta <- rbind(allMetaCor_pooled_bc, allMetaCor_dropone) %>%
  filter(metabType!="fapro") %>%
  arrange(-R2) %>% mutate(type=factor(type,levels=unique(type))) %>%
  select(R2, pval, type, metabType) %>%
  distinct() %>%
  ggplot() + geom_bar(aes(x=type, y=R2, fill=type), stat="identity", position="dodge") +
  facet_wrap(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("Dropped")
gg_R2drop_meta
ggsave("06a_metabolic_process/plot_pooledMRMcorrelations_droponeR2.png"
       ,gg_R2drop_meta)

gg_coefChange_meta <- rbind(allMetaCor_pooled_bc, allMetaCor_dropone) %>%
  filter(metabType!="fapro") %>%
  arrange(-R2) %>% mutate(type=factor(type,levels=unique(type))) %>%
  # select(R2, pval, type, metabType) %>%
  ggplot() + geom_point(aes(x=predictor, y=coef, col=type)) +
  facet_wrap(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("Dropped")
gg_coefChange_meta
ggsave("06a_metabolic_process/plot_pooledMRMcorrelations_coefChange.png"
       ,gg_coefChange_meta)


### Fapro ####
if ( !file.exists("06a_metabolic_process/allPCCorr_fapro.RData") ) {
  pointPCA_fapro <- as.data.frame(pca_fapro$x)
  allPCCorr_fapro <- data.frame()
  allPCCorr_fapro_pooled <- data.frame()
  for ( pc in 1:5 ) {
    # pc=1
    print(c("PC",pc))
    PCtemp <- as.matrix(pointPCA_fapro[,pc])
    rownames(PCtemp) <- rownames(pointPCA_fapro)
    PCtemp_dist <- as.matrix(dist(PCtemp))
    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist") ){
      # p="phyloDist"
      print(p)
      tempDist <- get(p)
      keepSites <- intersect(colnames(PCtemp_dist), colnames(tempDist))
      pcDist <- as.dist(PCtemp_dist[keepSites,keepSites])
      predDist <- as.dist(tempDist[keepSites,keepSites])
      mrm <- MRM(pcDist ~ predDist, mrank = TRUE, nperm = 1000)
      newdat <- data.frame(as.data.frame(mrm$coef) %>% rownames_to_column(var="predictor") %>%
                             filter(predictor!="Int") %>% 
                             mutate(predictor=p) %>%
                             rename(coef="pcDist", pval_coef=pval),t(mrm$r.squared)) %>%
        mutate(pc=pc)
      allPCCorr_fapro <- rbind(allPCCorr_fapro, newdat)
    }
    # Try all together
    listsites <- c()
    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
      tempDist <- as.matrix(get(p))
      listsites <- c(listsites,intersect(colnames(PCtemp_dist), colnames(tempDist)))
    }
    allKeepSites <- names(which(table(listsites)==5))

    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
      tempDist <- as.matrix(get(p))
      assign(paste0(p,"_filt"), as.dist(tempDist[allKeepSites,allKeepSites]))
    }
    PCtemp_dist_filt <- as.dist(PCtemp_dist[allKeepSites, allKeepSites])
    frml <- paste0("PCtemp_dist_filt ~ ", paste0(paste0(c("phyloDist","geoDist","ecoDist","climDist","timeDist"), "_filt"),collapse=" + "))
    mrm_temp <- MRM(formula(frml), mrank = TRUE)
    corrresults <- as.data.frame(mrm_temp$coef) %>% rownames_to_column(var="predictor") %>%
      rename(Coef="PCtemp_dist_filt", pval_coef=pval) %>% cbind(as.data.frame(t(mrm_temp$r.squared))) %>%
      mutate(type="fapro", pc=pc) %>% filter(predictor!="Int")
    allPCCorr_fapro_pooled <- rbind(corrresults, allPCCorr_fapro_pooled)
  }
  save(allPCCorr_fapro, file="06a_metabolic_process/allPCCorr_fapro.RData")
  save(allPCCorr_fapro_pooled, file="06a_metabolic_process/allPCCorr_fapro_pooled.RData")
  
} else {
  load("06a_metabolic_process/allPCCorr_fapro.RData")
  load("06a_metabolic_process/allPCCorr_fapro_pooled.RData")
}


gg_fapro <- allPCCorr_fapro %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=coef, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_fapro
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_fapro.png",
       gg_fapro)

gg_fapro_R2 <- allPCCorr_fapro %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=R2, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_fapro_R2
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_fapro_R2.png",
       gg_fapro_R2)

gg_fapropooled <- allPCCorr_fapro_pooled %>%
  filter(predictor!="Int") %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef, n=5, method="BH")) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=Coef, fill=predictor, col=padj<0.05), stat="identity", position="dodge", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_fapropooled
ggsave(filename = "06a_metabolic_process/Plot_allPCCorr_pooled_fapro.png"
       ,gg_fapropooled)


### ec ####
if ( ! file.exists("06a_metabolic_process/allPCCorr_ec.RData") ) {
  pointPCA_ec <- as.data.frame(pca_ec$x)
  
  allPCCorr_ec <- data.frame()
  allPCCorr_ec_pooled <- data.frame()
  for ( pc in 1:5 ) {
    print(c("PC",pc))
    PCtemp <- as.matrix(pointPCA_ec[,pc])
    rownames(PCtemp) <- rownames(pointPCA_ec)
    PCtemp_dist <- as.matrix(dist(PCtemp))
    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist") ){
      print(p)
      tempDist <- get(p)
      keepSites <- intersect(colnames(PCtemp_dist), colnames(tempDist))
      pcDist <- as.dist(PCtemp_dist[keepSites,keepSites])
      predDist <- as.dist(tempDist[keepSites,keepSites])
      mrm <- MRM(pcDist ~ predDist, mrank = TRUE, nperm = 1000)
      newdat <- data.frame(as.data.frame(mrm$coef) %>% rownames_to_column(var="predictor") %>%
                             filter(predictor!="Int") %>% 
                             mutate(predictor=p) %>%
                             rename(coef="pcDist", pval_coef=pval),t(mrm$r.squared)) %>%
        mutate(pc=pc)
      allPCCorr_ec <- rbind(allPCCorr_ec, newdat)
    }
    # Try all together
    listsites <- c()
    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
      tempDist <- as.matrix(get(p))
      listsites <- c(listsites,intersect(colnames(PCtemp_dist), colnames(tempDist)))
    }
    allKeepSites <- names(which(table(listsites)==5))

    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
      tempDist <- as.matrix(get(p))
      assign(paste0(p,"_filt"), as.dist(tempDist[allKeepSites,allKeepSites]))
    }
    PCtemp_dist_filt <- as.dist(PCtemp_dist[allKeepSites, allKeepSites])
    frml <- paste0("PCtemp_dist_filt ~ ", paste0(paste0(c("phyloDist","geoDist","ecoDist","climDist","timeDist"), "_filt"),collapse=" + "))
    mrm_temp <- MRM(formula(frml), mrank = TRUE)
    corrresults <- as.data.frame(mrm_temp$coef) %>% rownames_to_column(var="predictor") %>%
      rename(Coef="PCtemp_dist_filt", pval_coef=pval) %>% cbind(as.data.frame(t(mrm_temp$r.squared))) %>%
      mutate(type="ec", pc=pc)
    allPCCorr_ec_pooled <- rbind(corrresults, allPCCorr_ec_pooled)
  }
  # colnames(allPCCorr_ec) <-  c("mantelr","pval1","pval2","pval3","llim.2.5%","ulim.97.5%","type","pc" )
  save(allPCCorr_ec, file="06a_metabolic_process/allPCCorr_ec.RData")
  save(allPCCorr_ec_pooled, file="06a_metabolic_process/allPCCorr_ec_pooled.RData")
} else {
  load("06a_metabolic_process/allPCCorr_ec.RData")
  load("06a_metabolic_process/allPCCorr_ec_pooled.RData")
}

gg_ec <- allPCCorr_ec %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=coef, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_ec
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_ec.png",
       gg_ec)

gg_ec_R2 <- allPCCorr_ec %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=R2, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_ec_R2
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_ec_R2.png",
       gg_ec_R2)

gg_ecpooled <- allPCCorr_ec_pooled %>%
  filter(predictor!="Int") %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef, n=5, method="BH")) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=Coef, fill=predictor, col=padj<0.05), stat="identity", position="dodge", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_ecpooled
ggsave(filename = "06a_metabolic_process/Plot_allPCCorr_pooled_ec.png"
       ,gg_ecpooled)

### ko ####
if ( !file.exists("06a_metabolic_process/allPCCorr_ko.RData") ) {
  pointPCA_ko <- as.data.frame(pca_ko$x)
  
  allPCCorr_ko <- data.frame()
  allPCCorr_ko_pooled <- data.frame()
  for ( pc in 1:5 ) {
    print(c("PC",pc))
    PCtemp <- as.matrix(pointPCA_ko[,pc])
    rownames(PCtemp) <- rownames(pointPCA_ko)
    PCtemp_dist <- as.matrix(dist(PCtemp))
    for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist") ){
      print(p)
      tempDist <- get(p)
      keepSites <- intersect(colnames(PCtemp_dist), colnames(tempDist))
      pcDist <- as.dist(PCtemp_dist[keepSites,keepSites])
      predDist <- as.dist(tempDist[keepSites,keepSites])
      mrm <- MRM(pcDist ~ predDist, mrank = TRUE, nperm = 1000)
      newdat <- data.frame(as.data.frame(mrm$coef) %>% rownames_to_column(var="predictor") %>%
                             filter(predictor!="Int") %>% 
                             mutate(predictor=p) %>%
                             rename(coef="pcDist", pval_coef=pval),t(mrm$r.squared)) %>%
        mutate(pc=pc)
      allPCCorr_ko <- rbind(allPCCorr_ko, newdat)
    }
    # # Try all together
    # listsites <- c()
    # for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
    #   tempDist <- as.matrix(get(p))
    #   listsites <- c(listsites,intersect(colnames(PCtemp_dist), colnames(tempDist)))
    # }
    # allKeepSites <- names(which(table(listsites)==5))
    # 
    # for ( p in c("phyloDist","geoDist","ecoDist","climDist","timeDist")  ) {
    #   tempDist <- as.matrix(get(p))
    #   assign(paste0(p,"_filt"), as.dist(tempDist[allKeepSites,allKeepSites]))
    # }
    # PCtemp_dist_filt <- as.dist(PCtemp_dist[allKeepSites, allKeepSites])
    # frml <- paste0("PCtemp_dist_filt ~ ", paste0(paste0(c("phyloDist","geoDist","ecoDist","climDist","timeDist"), "_filt"),collapse=" + "))
    # mrm_temp <- MRM(formula(frml), mrank = TRUE)
    # corrresults <- as.data.frame(mrm_temp$coef) %>% rownames_to_column(var="predictor") %>%
    #   rename(Coef="PCtemp_dist_filt", pval_coef=pval) %>% cbind(as.data.frame(t(mrm_temp$r.squared))) %>%
    #   mutate(type="ko", pc=pc)
    # allPCCorr_ko_pooled <- rbind(corrresults, allPCCorr_ko_pooled)
  }
  # colnames(allPCCorr_ko) <-  c("mantelr","pval1","pval2","pval3","llim.2.5%","ulim.97.5%","type","pc" )
  save(allPCCorr_ko, file="06a_metabolic_process/allPCCorr_ko.RData")
  # save(allPCCorr_ko_pooled, file="06a_metabolic_process/allPCCorr_ko_pooled.RData")
} else {
  load("06a_metabolic_process/allPCCorr_ko.RData")
  load("06a_metabolic_process/allPCCorr_ko_pooled.RData")
}

gg_ko <- allPCCorr_ko %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=coef, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_ko
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_ko.png",
       gg_ko)

gg_ko_R2 <- allPCCorr_ko %>%
  rowwise() %>%
  mutate(padj_coef=p.adjust(pval_coef, method="BH", n = 5)) %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=R2, fill=predictor, col=padj_coef<0.05), position="dodge", stat="identity", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_ko_R2
ggsave(filename="06a_metabolic_process/Plot_allPCCorr_ko_R2.png",
       gg_ko_R2)


gg_kopooled <- allPCCorr_ko_pooled %>%
  filter(predictor!="Int") %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef, n=5, method="BH")) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x=factor(pc), y=Coef, fill=predictor, col=padj<0.05), stat="identity", position="dodge", width = 0.5) +
  scale_color_manual(values=c("white","black"))
gg_kopooled
ggsave(filename = "06a_metabolic_process/Plot_allPCCorr_pooled_ko.png"
       ,gg_kopooled)

#### total loadings of PCs and enviro #####
# allPCCorr_ko  %>%
  # ggplot() + geom_bar(aes(x=factor(pc), y=mantelr, fill=type), position="stack", stat="identity") 

write.table(allPCCorr_ec, file="06a_metabolic_process/table_pc_corr_env_ec.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(allPCCorr_ko, file="06a_metabolic_process/table_pc_corr_env_ko.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(allPCCorr_fapro, file="06a_metabolic_process/table_pc_corr_env_fapro.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(allPCCorr_ec_pooled, file="06a_metabolic_process/table_pc_corr_env_ec.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(allPCCorr_ko_pooled, file="06a_metabolic_process/table_pc_corr_env_ko.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(allPCCorr_fapro_pooled, file="06a_metabolic_process/table_pc_corr_env_fapro.txt", quote=FALSE, row.names=FALSE, sep="\t")


## Variance explained
EC_PC_variances <- (pca_ec$sdev^2)/sum(pca_ec$sdev^2)
KO_PC_variances <- (pca_ko$sdev^2)/sum(pca_ko$sdev^2)
FAPRO_PC_variances <- (pca_fapro$sdev^2)/sum(pca_fapro$sdev^2)

sink("06a_metabolic_process/Variation_by_PC.txt")
print("fapro first five PC")
FAPRO_PC_variances[c(1:5)]
print("EC first five PC")
EC_PC_variances[c(1:5)]
print("EC first five PC")
KO_PC_variances[c(1:5)]
print("sum of first five PCs EC")
sum(EC_PC_variances[c(1:5)])
print("sum of first five PCs KO")
sum(KO_PC_variances[c(1:5)])
sink()


### Custom fapro plot #####
topFunc <- names(sort(pca_fapro$center, decreasing = T)[1:8])
gg_fapro_firsfivePC <- as.data.frame(pca_fapro$rotation) %>% rownames_to_column(var="Function") %>%
  select(Function, PC1, PC2, PC3, PC4, PC5) %>% 
  filter(Function %in% topFunc) %>% pivot_longer(-Function, names_to="Principal component", values_to="Correlation with PC") %>%
  ggplot() +geom_bar(aes(x=Function, y=`Correlation with PC`, fill=Function), stat="identity", position="dodge", show.legend = FALSE) +
  facet_grid(`Principal component`~.) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
gg_fapro_firsfivePC
ggsave(filename = "06a_metabolic_process/fapro_first5_PC.png", height=7, width=3, gg_fapro_firsfivePC)
