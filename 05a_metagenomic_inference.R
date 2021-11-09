#!bin/bash Rscript

### Metagenomic analysis ####
library(vegan)
library(ggfortify)
library(ecodist)
library(tidyverse)

## Load ##
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt") 
load("04d_host_and_enviro_matrices/RData/climDist.RData")
climDist <- as.matrix(climDist)
load("04d_host_and_enviro_matrices/RData/ecoDist.RData")
ecoDist <- as.matrix(ecoDist)
load("04d_host_and_enviro_matrices/RData/geoDist.RData")
geoDist <- as.matrix(geoDist)
load("04d_host_and_enviro_matrices/RData/phyloDist.RData")
phyloDist <- as.matrix(phyloDist)
load("04d_host_and_enviro_matrices/RData/timeDist.RData")
timeDist <- as.matrix(timeDist)
load("04d_host_and_enviro_matrices/RData/studyDist.RData")
studyDist <- as.matrix(studyDist)

fapro <- read.delim("04c_FAPROTAX/func_table.tsv") 
fapro <- fapro[, c("X.group", (meta[meta$sampleid %in% colnames(fapro),"sampleid"]))]
rnfapro <- fapro$X.group
fapro <- fapro[,-which(colnames(fapro)=="X.group")]
rownames(fapro)=rnfapro

ec <- read.delim("02b_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv")
ec <- ec[,c(c("function.","description", (meta[meta$sampleid %in% colnames(ec),"sampleid"])) )]
ec_ref <- ec[,c("function.","description")]
colnames(ec_ref) <- c("MetabolicGroup","description")
ec <- ec[,-which(colnames(ec)%in%c("function.","description"))]
rownames(ec) <- ec_ref$MetabolicGroup

ko <- read.delim("02b_picrust/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv")
ko <- ko[,c(c("function.","description", (meta[meta$sampleid %in% colnames(ko),"sampleid"])) )]
ko_ref <- ko[,c("function.","description")]
colnames(ko_ref) <- c("MetabolicGroup","description")
ko <- ko[,-which(colnames(ko)%in%c("function.","description"))]
rownames(ko) <- ko_ref$MetabolicGroup

#### 
dir.create("05a_metagenomic_inference")
dir.create("05a_metagenomic_inference/lefse_data")
dir.create("05a_metagenomic_inference/fapro_plots")
dir.create("05a_metagenomic_inference/ec_plots")
dir.create("05a_metagenomic_inference/ko_plots")

source("code/metagenomic_inference_functions.R")
# First, let's look at the functions that define samples
# Most common in dataset are:
### List most abundant functions



#### Run and save summaries OR load summaries ####

if ( ! file.exists("05a_metagenomic_inference/koMetab.RData")) {
  ### FAPRO
  fapro_filt <- fapro[rowSums(fapro)>0,colSums(fapro)>0]
  # Most common function
  fapro_observedInstances <- data.frame(ObservedInstances=sort(rowSums(fapro_filt))) %>% as.data.frame() %>%
    rownames_to_column(var="metab") %>%
    mutate(ObservedInstances_rel=ObservedInstances/sum(ObservedInstances)) %>%
    arrange(-ObservedInstances)
  write.table(fapro_observedInstances, file = "05a_metagenomic_inference/fapro_observed_instances.txt", quote=FALSE, row.names=FALSE, sep="\t")
  

  # If already ran, re-load
  if ( file.exists("05a_metagenomic_inference/pca_fapro_filt_relabund.RData")) {
    load("05a_metagenomic_inference/pca_fapro_filt_relabund.RData")
    importPCA <- pca_fapro
  } else {
    importPCA <- NULL
  }
  fapro_summary <- summarise_functional_table(tab=fapro_filt, metab_single_cutoff = 0.005
                                              , importPCA=importPCA)
  # If haven't saved pca, save it now
  if (! file.exists("05a_metagenomic_inference/pca_fapro_filt_relabund.RData")) {
    pca_fapro <- fapro_summary$pca
    save(pca_fapro, file="05a_metagenomic_inference/pca_fapro_filt_relabund.RData")}
  
  ggsave(filename="05a_metagenomic_inference/fapro_plots/pca.png", height=6, width=8
         ,fapro_summary$Plot_pca )
  ggsave(filename="05a_metagenomic_inference/fapro_plots/PC_by_metab.png", height=6, width=10
         ,fapro_summary$Plot_PC_by_metab )
  ggsave(filename="05a_metagenomic_inference/fapro_plots/PC_var.png", height=3, width=4
         ,fapro_summary$Plot_PC_Variance_top )
  ggsave(filename="05a_metagenomic_inference/fapro_plots/metab_top.png", height=8, width=8
         ,fapro_summary$Plot_metabsummary_top )
  fapro_pc_var <- data.frame(PC_variation=(fapro_summary$pca$sdev^2)/sum(fapro_summary$pca$sdev^2), PC_cumsum=cumsum((fapro_summary$pca$sdev^2)/sum(fapro_summary$pca$sdev^2)))
  write.table(fapro_pc_var, file="05a_metagenomic_inference/fapro_PC_variationexplained.txt", quote=FALSE, row.names = FALSE, sep="\t")
  
  ## EC
  ec_filt <- ec[rowSums(ec)>0, colSums(ec)>0]
  ec_top <- apply(ec_filt, 2, function(x) x/sum(x)) %>%
    rowMeans() %>% data.frame() %>% rownames_to_column(var="MetabolicGroup") %>%
    rename(meanRelativeAbund=".") %>% left_join(ec_ref) %>%
    arrange(-meanRelativeAbund)
  write.table(ec_top, file="05a_metagenomic_inference/ec_topmetab.txt", quote=FALSE, row.names = FALSE, sep="\t")
  
  if ( file.exists("05a_metagenomic_inference/pca_ec_filt_relabund.RData")) {
    load("05a_metagenomic_inference/pca_ec_filt_relabund.RData")
    importPCA <- pca_ec
  } else {
    importPCA <- NULL
  }
  ec_summary <- summarise_functional_table(tab=ec_filt, metab_single_cutoff = 0.05
                                           , ref_labs = ec_ref
                                           , ref_col_names = colnames(ec_ref)
                                           , importPCA=importPCA)
  # If haven't saved pca, save it now
  if (! file.exists("05a_metagenomic_inference/pca_ec_filt_relabund.RData")) {
    pca_ec <- ec_summary$pca
    save(pca_ec, file="05a_metagenomic_inference/pca_ec_filt_relabund.RData")}
  
  
  ggsave(filename="05a_metagenomic_inference/ec_plots/pca.png", height=6, width=8
         ,ec_summary$Plot_pca )
  ggsave(filename="05a_metagenomic_inference/ec_plots/PC_by_metab.png", height=6, width=10
         ,ec_summary$Plot_PC_by_metab )
  ggsave(filename="05a_metagenomic_inference/ec_plots/PC_var.png", height=3, width=4
         ,ec_summary$Plot_PC_Variance_top )
  ggsave(filename="05a_metagenomic_inference/ec_plots/metab_top.png", height=8, width=8
         ,ec_summary$Plot_metabsummary_top )
  ec_pc_var <- data.frame(PC_variation=(ec_summary$pca$sdev^2)/sum(ec_summary$pca$sdev^2), PC_cumsum=cumsum((ec_summary$pca$sdev^2)/sum(ec_summary$pca$sdev^2)))
  write.table(ec_pc_var[1:100,], file="05a_metagenomic_inference/ec_PC_variationexplained.txt", quote=FALSE, row.names = FALSE, sep="\t")

  
  #### KO 
  ko_filt <- ko[rowSums(ko)>0, colSums(ko)>0]
  ko_top <- apply(ko_filt, 2, function(x) x/sum(x)) %>%
    rowMeans() %>% data.frame() %>% rownames_to_column(var="MetabolicGroup") %>%
    rename(meanRelativeAbund=".") %>% left_join(ko_ref) %>%
    arrange(-meanRelativeAbund)
  write.table(ko_top[1:100,], file="05a_metagenomic_inference/ec_topmetab.txt", quote=FALSE, row.names = FALSE, sep="\t")
  if ( file.exists("05a_metagenomic_inference/pca_ko_filt_relabund.RData")) {
    load("05a_metagenomic_inference/pca_ko_filt_relabund.RData")
    importPCA <- pca_ko
  } else {
    importPCA <- NULL
  }
  ko_summary <- summarise_functional_table(tab=ko_filt, metab_single_cutoff = 0.10
                                           , ref_labs = ko_ref
                                           , ref_col_names = colnames(ko_ref)
                                           , importPCA=importPCA)
  # If haven't saved pca, save it now
  if (! file.exists("05a_metagenomic_inference/pca_ko_filt_relabund.RData")) {
    pca_ko <- ko_summary$pca
    save(pca_ko, file="05a_metagenomic_inference/pca_ko_filt_relabund.RData")}
  
  ggsave(filename="05a_metagenomic_inference/ko_plots/pca.png", height=6, width=8
         ,ko_summary$Plot_pca )
  ggsave(filename="05a_metagenomic_inference/ko_plots/PC_by_metab.png", height=6, width=10
         ,ko_summary$Plot_PC_by_metab )
  ggsave(filename="05a_metagenomic_inference/ko_plots/PC_var.png", height=3, width=4
         ,ko_summary$Plot_PC_Variance_top )
  ggsave(filename="05a_metagenomic_inference/ko_plots/metab_top.png", height=8, width=8
         ,ko_summary$Plot_metabsummary_top )
  ko_pc_var <- data.frame(PC_variation=(ko_summary$pca$sdev^2)/sum(ko_summary$pca$sdev^2), PC_cumsum=cumsum((ko_summary$pca$sdev^2)/sum(ko_summary$pca$sdev^2)))
  write.table(ko_pc_var[1:100,], file="05a_metagenomic_inference/ko_PC_variationexplained.txt", quote=FALSE, row.names = FALSE, sep="\t")
  
  ### SAVE ####
  
  faproMetab <- fapro_summary$top_metab
  save(faproMetab, file="05a_metagenomic_inference/faproMetab.RData")
  
  ecMetab <- ec_summary$top_metab
  save(ecMetab, file="05a_metagenomic_inference/ecMetab.RData")
  
  koMetab <- ko_summary$top_metab
  save(koMetab, file="05a_metagenomic_inference/koMetab.RData")
} else {
  load("05a_metagenomic_inference/faproMetab.RData")
  load("05a_metagenomic_inference/ecMetab.RData")
  load("05a_metagenomic_inference/koMetab.RData")
}

### Correlate with dms
if ( ! file.exists("05a_metagenomic_inference/allMetaCor_bc.RData") ) {
  # allMetaCor <- data.frame()
  # allMetaCor_pooled <- data.frame()
  allMetaCor_bc <- data.frame()
  allMetaCor_pooled_bc <- data.frame()
  for ( mt in c("fapro","ec","ko") ) {
    # mt="fapro"
    print(paste0("doing ",mt))
    toKeepMetab <- get(paste0(mt,"Metab"))
    tempFilt <- get(paste0(mt))[toKeepMetab,]
    tempFilt <- tempFilt[,colSums(tempFilt)>0]
    tempFilt <- apply(tempFilt, 2, function(x) x/sum(x))
    # tempDist <- as.matrix(dist(t(tempFilt)))
    tempDist <- as.matrix(vegdist(t(tempFilt), method="bray"))
    for ( ed in c("phyloDist","climDist","geoDist","ecoDist","timeDist") ) {
      # ed="phyloDist"
      print(paste0("type ",ed))
      
      EnvDist <- get(ed)
      sharedSites <- intersect(colnames(tempFilt), colnames(EnvDist))
      tempDist_filt <- as.dist(tempDist[sharedSites,sharedSites])
      EnvDist <- as.dist(EnvDist[sharedSites,sharedSites])
      # Make dist
      testCorr <- MRM(tempDist_filt ~ EnvDist, mrank = TRUE)
      newdat <- data.frame(as.data.frame(testCorr$coef) %>% rownames_to_column(var="predictor") %>%
        filter(predictor!="Int") %>% rename(coef=tempDist_filt, pval_coef=pval),t(testCorr$r.squared)) %>%
        mutate(type=ed, metabType=mt)
      # allMetaCor <- rbind(newdat, allMetaCor)
      allMetaCor_bc <- rbind(newdat, allMetaCor_bc)
    }
    
    ## Do them pooled
    allsites <- c()
    for ( ed in c("phyloDist","climDist","geoDist","ecoDist","timeDist")) {
      tempPred <- as.matrix(get(paste0(ed)))
      allsites <- c(allsites, intersect(colnames(tempPred), colnames(tempFilt)))
    }
    filtSites <- names(which(table(allsites)==5))
    for ( ed in c("phyloDist","climDist","geoDist","ecoDist","timeDist")) {
      tempPred <- as.matrix(get(paste0(ed)))
      assign(paste0(ed,"_filt"), as.dist(tempPred[filtSites,filtSites]))
    }    
    tempDist_filt <- as.dist(tempDist[filtSites,filtSites])
    frml <- paste0("tempDist_filt ~ ",paste0(paste0(c("phyloDist","climDist","geoDist","ecoDist","timeDist"), "_filt"), collapse=" + "))
    mrm_pooled <- MRM(formula(frml), mrank = T, nperm = 1000)
    newdat_pooled <- data.frame(as.data.frame(mrm_pooled$coef) %>% rownames_to_column(var="predictor") %>%
                           filter(predictor!="Int") %>% rename(coef=tempDist_filt, pval_coef=pval),t(mrm_pooled$r.squared)) %>%
      mutate(type="pooled", metabType=mt)
    # allMetaCor_pooled <- rbind(newdat_pooled, allMetaCor_pooled)
    allMetaCor_pooled_bc <- rbind(newdat_pooled, allMetaCor_pooled_bc)
  }
  # save(allMetaCor, file = "05a_metagenomic_inference/allMetaCor.RData")
  # save(allMetaCor_pooled, file = "05a_metagenomic_inference/allMetaCor_pooled.RData")
  save(allMetaCor_bc, file = "05a_metagenomic_inference/allMetaCor_bc.RData")
  save(allMetaCor_pooled_bc, file = "05a_metagenomic_inference/allMetaCor_pooled_bc.RData")
  
} else {
  # load("05a_metagenomic_inference/allMetaCor.RData")
  # load("05a_metagenomic_inference/allMetaCor_pooled.RData")
  load("05a_metagenomic_inference/allMetaCor_bc.RData")
  load("05a_metagenomic_inference/allMetaCor_pooled_bc.RData")
}


### DO DROP ONE ###
listPred <- c("phyloDist","climDist","geoDist","ecoDist","timeDist")
if ( ! file.exists("05a_metagenomic_inference/allMetaCor_dropone.RData") ) {
  allMetaCor_dropone <- data.frame()
  for ( mt in c("ec","ko") ) {
    print(paste0("doing ",mt))
    toKeepMetab <- get(paste0(mt,"Metab"))
    tempFilt <- get(paste0(mt))[toKeepMetab,]
    tempFilt <- tempFilt[,colSums(tempFilt)>0]
    tempFilt <- apply(tempFilt, 2, function(x) x/sum(x))
    # tempDist <- as.matrix(dist(t(tempFilt)))
    tempDist <- as.matrix(vegdist(t(tempFilt), method="bray"))
    ## Do them pooled
    allsites <- c()
    for ( ed in listPred) {
      tempPred <- as.matrix(get(paste0(ed)))
      allsites <- c(allsites, intersect(colnames(tempPred), colnames(tempFilt)))
    }
    filtSites <- names(which(table(allsites)==5))
    for ( ed in listPred) {
      tempPred <- as.matrix(get(paste0(ed)))
      assign(paste0(ed,"_filt"), as.dist(tempPred[filtSites,filtSites]))
    }    
    tempDist_filt <- as.dist(tempDist[filtSites,filtSites])
    # DROP ONE
    for ( d in listPred ) {
      newList <- listPred[listPred!=d]
      frml <- paste0("tempDist_filt ~ ",paste0(paste0(newList, "_filt"), collapse=" + "))
      mrm_pooled <- MRM(formula(frml), mrank = T, nperm = 1000)
      newdat_pooled <- data.frame(as.data.frame(mrm_pooled$coef) %>% rownames_to_column(var="predictor") %>%
                                    filter(predictor!="Int") %>% rename(coef=tempDist_filt, pval_coef=pval),t(mrm_pooled$r.squared)) %>%
        mutate(type=paste0("drop_",d), metabType=mt)
      # allMetaCor_pooled <- rbind(newdat_pooled, allMetaCor_pooled)
      allMetaCor_dropone <- rbind(newdat_pooled, allMetaCor_dropone)
    }
    
  }
  save(allMetaCor_dropone, file = "05a_metagenomic_inference/allMetaCor_dropone.RData")
} else {
  load("05a_metagenomic_inference/allMetaCor_dropone.RData")
}

# 
# gg_metagCor <- allMetaCor %>%
#   ggplot() + geom_point(aes(x=type, y=R2)) +
#   facet_wrap(.~metabType) +
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
# 
# ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_basic.png"
#        ,gg_metagCor)


gg_metagCor <- allMetaCor_bc %>%
  ggplot() + geom_point(aes(x=type, y=R2)) +
  facet_wrap(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_basic_bc.png"
       ,gg_metagCor)


# 
# gg_metagCor_pooled <- allMetaCor_pooled %>%
#   ggplot() + geom_point(aes(x=predictor, y=coef)) +
#   facet_wrap(.~metabType) +
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
# 
# ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_basic_pooled.png"
#        ,gg_metagCor_pooled)


gg_metagCor_pooled <- allMetaCor_pooled_bc %>%
  ggplot() + geom_point(aes(x=predictor, y=coef)) +
  facet_wrap(.~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_basic_pooled_bc.png"
       ,gg_metagCor_pooled)


gg_metagCor_dropone_coef <- allMetaCor_dropone %>% 
  ggplot() + geom_bar(aes(x=predictor, y=coef), stat='identity') +
  facet_grid(type~metabType) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
gg_metagCor_dropone_coef
ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_dropone.png"
       ,gg_metagCor_dropone_coef)

rbind(allMetaCor_pooled_bc, allMetaCor_dropone) %>% 
  filter(metabType!="fapro") %>%
  select(R2, type, metabType) %>% distinct() %>%
  ggplot() + geom_jitter(aes(x=metabType, y=R2, col=type), stat='identity') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) 
allMetaCor_dropone_R2
ggsave(filename = "05a_metagenomic_inference/metagonme_correlation_dropone.png"
       ,allMetaCor_dropone_R2)
