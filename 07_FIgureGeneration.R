#!bin/bash
# library(castor)
# library(caper)
library(ape)
library(ggfortify)
library(tidyverse)
library(cowplot)

###### FIGURE GENERATION ###########

dir.create("07_FigureGeneration")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt") 
amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
amphtree <- keep.tip(amphtree, unique(meta$Host))


###### Predictor data ######
#### %--------Bioclimate ####
env_fact <- c("biom_AnnualMeanTemp", "biom_MeanDiurnalTempRange", "biom_TempSeasonality", "biom_TempAnnualRange", "biom_AnnualPrecip", "biom_PrecipSeasonality"  )
meta_climate <- meta %>% select(Longitude, Latitude, one_of(env_fact))
# meta_climate <- meta %>% select(Longitude, Latitude, starts_with("biom_")) %>% distinct() %>% drop_na() %>%
# rename_at(vars(starts_with("biom_")), ~gsub("biom_","",.))
  
# meta_climate[,3]-mean(meta_climate[,3])
meta_climate_st <- apply(meta_climate[,-c(1,2)],2, function(x) (x-mean(x))/sd(x))
pca_env <- prcomp(meta_climate_st)
gg_climatepcoa <-  autoplot(pca_env, data=meta_climate_st
                     , loadings=TRUE
                     , loadings.label=TRUE) 
# arrow_ends <- layer_data(gg_temp,2)
# rownames(arrow_ends) <- names(pca_env$center)
# textClim <- arrow_ends[c("AnnualPrecip","TempSeasonality"),] %>%
#   rownames_to_column(var="Biom") %>% as.data.frame()
# 
# gg_climatepcoa <- gg_temp+
#   geom_text(textClim, mapping=aes(x=xend, y=yend, label=Biom), nudge_y = 0.03)
gg_climatepcoa +
  xlim(c(-0.05, 0.05))
ggsave(filename="07_FigureGeneration/predictor_PCA_bioms.png", width=6.5, height=6
       ,gg_climatepcoa+
         xlim(c(-0.05, 0.05)))


#### %--------Geographic location #####
# Plot map
latlong <- meta[,c("Latitude","Longitude")] %>% distinct()
# Get number of samples
nSampDat <- meta %>% group_by(Latitude, Longitude) %>% summarise(nSamp=n()) %>% ungroup() %>%
  mutate(size=(nSamp+20)*10/max(nSamp+20)) 
nSampDat <- left_join(latlong, nSampDat)

png("07_FigureGeneration/map_samples_sizeSamp.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="white", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg="yellow", pch=21, cex=nSampDat$size)
dev.off()


### %--------Amphibian Tree ####
strHab <- meta$HabitatClass[match(amphtree$tip.label, meta$Host)]
habcols <- data.frame(Habitat=c("Aquatic","Semi-aquatic","Arboreal","Terrestrial"), habCol=c("blue","purple","darkgreen","darkred"))
strHabCol <- left_join(data.frame(Habitat=strHab),habcols) %>% pull(habCol)
png("07_FigureGeneration/amphibians_with_eco.png", height=1200, width = 800)
plot(amphtree, tip.color=strHabCol, cex = 1.8)
legend("topleft", c(habcols$Habitat), col= habcols$habCol, pch=19, cex=1.6)
dev.off()


### %--------Time of samples ####
meta_calendar <- meta %>% select(sampleid, Host, Location, Month)%>%
  # filter(Month%in% c(1,2,3,4,5,6,7,8,9,10,11,)) %>%
  group_by(Host, Month) %>% summarise(NumberSamples=n()) %>% ungroup() %>%
  arrange(Host) %>% 
  mutate(HostNumber=group_indices(., Host)/65) %>%
  mutate(theta=pi/2-2*pi*(Month-1)/12, H=0.5+HostNumber) %>%
  mutate(O_Y=H*sin(theta), A_X=H*cos(theta)) 
calendar <- data.frame(Month=seq(1,12), MonthText=c("J","F","Mr","Ap","My","Jn","Jl","Au","S","O","N","D")) %>% mutate(theta=pi/2-2*pi*(Month-1)/12, H=0.4) %>%
  mutate(O_Y=H*sin(theta), A_X=H*cos(theta), yend=H*0.8*sin(theta), xend=H*0.8*cos(theta))

gg_calendar <- ggplot(meta_calendar) + geom_point(aes(x=A_X, y=O_Y, cex=NumberSamples, col=factor(Host)))+
  geom_text(data=calendar, aes(x=A_X, y=O_Y, label=MonthText))+
  geom_segment(data=calendar, aes(x=0,xend=xend, y=0, yend=yend)) +
  xlim(-1.5,1.5) + ylim(-1.5,1.5) +
  geom_point(data.frame(x=0,y=0), mapping=aes(x=x, y=y)) +
  theme(axis.text=element_blank(), axis.line = element_blank(), axis.ticks = element_blank()
        , rect = element_rect(), axis.title = element_blank(), panel.background = element_rect(fill="white")
        , panel.grid = element_blank(), panel.border = element_rect(colour='black', fill=NA)) +
  xlab("")+ ylab("") +
  guides(colour="none") + labs(cex="Number of\nSamples")
gg_calendar
ggsave(filename = "07_FigureGeneration/calendar_samples.png", height=4, width=6,
       gg_calendar)


######### Mapping colors ########
mapPred <- data.frame(predictor=c("pooled","phyloDist_filt","geoDist_filt","climDist_filt","timeDist_filt","studyDist_filt","ecoDist_filt")
                      , Predictor=c("All predictors together","Phylogeny\n(Amphibian phylo. distance)","Geography\n(Lat/Long distance)","Climate\n(Difference in climate)","Time\n(Diff in sampling month)","Study\n(Effect of study)","Ecology\n(Arb, Terr, Aq, Semi-Aq)")
                      , drop = c("fulldataset", "nophylo", "nogeo","noclim","notime","nostudy","noeco")
                      , Model= c("Full dataset","Phylogeny omitted","Geography omitted","Climate omitted","Time omitted","Study omitted","Ecology omitted")
                      , PredCol = c("black","blue","darkgreen","red","purple","pink","lightgreen")
                      , Simple=c("All predictors","Phylogeny","Geography","Climate","Time","Study","Ecology"))
mapPred <- mapPred[match(c("pooled","geoDist_filt","climDist_filt","ecoDist_filt","timeDist_filt","phyloDist_filt","studyDist_filt"),mapPred$predictor),]
predColors <- mapPred$PredCol
names(predColors) <- mapPred$Predictor
dropColors <- mapPred$PredCol
names(dropColors) <- mapPred$Model


######### BBDT plots #########
load("06a_BDTT_process/allDat_pooled_withdrop.RData")
load("06a_BDTT_process/allDat_bdtt_separate.RData")
just_pooled <- allDat_pooled_withdrop%>%filter(drop=="fulldataset") %>% select(R2,pval,slice)%>%distinct() %>%mutate(predictor="pooled")
allDat_bdtt_withpooled <- allDat_bdtt %>% 
  full_join(just_pooled) 

#### %--------Separate correlations R2 #####

# Map predictors

gg_BDTT_separateR2 <- allDat_bdtt_withpooled %>%
  full_join(mapPred) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval,n=6)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  ggplot(aes(x=slice, y=R2, col=Predictor)) + 
  geom_line() +
  geom_point(aes(pch=Significant), cex=2) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),legend.key.size = unit(0.75,"cm"), legend.text = element_text(size=8))+
  scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  scale_color_manual(values= predColors) +
  ylab(expression(Variance~explained~by~predictor~(R^2))) +
  xlab("Depth into phylogenetic tree\nat which tips are collapsed")
gg_BDTT_separateR2
ggsave(filename = "07_FigureGeneration/BDTT_separate_R2.png", height=5, width=7
       ,gg_BDTT_separateR2)

##### %--------Text R2 drops #####
maxR2bdtt <- allDat_pooled_withdrop %>% filter(drop=="fulldataset") %>%
  select(R2, slice) %>%
  rename(R2_total=R2) %>% distinct()

tab_bdttR2drop <- allDat_pooled_withdrop %>%
  select(R2, drop, slice) %>% distinct() %>%
  left_join(maxR2bdtt) %>%
  mutate(R2_drop=R2_total-R2) %>%
  select(slice, drop, R2_drop) %>%
  filter(slice<0.5) %>%
  pivot_wider(names_from=slice, values_from=R2_drop)
tab_bdttR2drop
write.table(tab_bdttR2drop, file="07_FigureGeneration/table_BDTT_R2drop.txt", quote=FALSE, row.names = FALSE, sep="\t")

##### %--------Text R2 0.05 #####
R20.05 <- allDat_bdtt %>% select(R2, predictor, slice) %>% 
  filter(slice==0.05) %>% distinct()
tab_bdttR20.5 <- allDat_pooled_withdrop %>%
  filter(drop=="fulldataset") %>%
  filter(slice==0.05) %>%
  select(R2, drop, slice) %>% distinct() %>%
  rename(predictor=drop) %>%
  full_join(R20.05) %>%
  arrange(-R2) 
tab_bdttR20.5
write.table(tab_bdttR20.5, file="07_FigureGeneration/table_BDTT_R20.05.txt", quote=FALSE, row.names = FALSE, sep="\t")


##### %--------One with total R2 of each drop-one out #####
gg_BDTT_drop <- allDat_pooled_withdrop %>%
  select(R2, pval, slice,drop) %>% distinct() %>%
  full_join(mapPred) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval,n=6)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  ggplot(aes(x=slice, y=R2, col=Model)) + 
  geom_point(aes(pch=Significant), cex=2, position=position_dodge(width=0.05), alpha=0.5) +
  geom_line(position=position_dodge(width=0.05), alpha=0.5) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_manual(values= dropColors) +
  scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  ylab(expression(Variance~explained~by~model~(R^2))) +
  xlab("Depth into phylogenetic tree\nat which tips are collapsed")
gg_BDTT_drop
ggsave(filename = "07_FigureGeneration/BDTT_droptest.png", height=4.5, width=7
       ,gg_BDTT_drop)

#### %-------- Combined R2 ####
gg_aligned_bdtt <- plot_grid(gg_BDTT_separateR2 + labs(title="a."), gg_BDTT_drop+labs(title="b."), ncol=1, align="v")
gg_aligned_bdtt
ggsave(filename = "07_FigureGeneration/BDTT_two_panel_figure.png", height=8, width=7,
       gg_aligned_bdtt)

##### %--------Coef of drop-one comparison #####
gg_BDTT_drop_coef <- allDat_pooled_withdrop %>%
  mutate(predmat = ifelse(predictor=="climDist_filt", "Climate predictor",
                          ifelse(predictor=="ecoDist_filt","Host ecology predictor",
                                 ifelse(predictor=="geoDist_filt", "Geographic distance predictor",
                                        ifelse(predictor=="phyloDist_filt","Host phylogeny predictor"
                                               , ifelse(predictor=="studyDist_filt","Study predictor",
                                                        ifelse(predictor=="timeDist_filt","Time predictor", NA))))))) %>%
  
  select(predmat, Coef, pval_coef, slice, drop) %>% 
  # select(R2, pval, slice,drop) %>% distinct() %>%
  full_join(mapPred %>% select(-predictor)) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef,n=6)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  arrange(Model) %>% mutate(predmat=factor(predmat, levels=unique(predmat))) %>%
  ggplot(aes(x=slice, y=abs(Coef), col=Model)) + 
  geom_point(aes(pch=Significant), cex=2, position=position_dodge(width=0.05)) +
  geom_line(position=position_dodge(width=0.05), alpha=0.5) +
  facet_wrap(.~predmat) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_manual(values= dropColors) +
  scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  ylab("Regression coefficient") +
  xlab("Depth into phylogenetic tree\nat which tips are collapsed")
  gg_BDTT_drop_coef
ggsave(filename = "07_FigureGeneration/BDTT_droptest_coeff.png", height=5, width=8
       ,gg_BDTT_drop_coef)

#### Depth of tip collapsed stats ####
load("06a_BDTT_process/allSlicesStats_bdtt.RData")
####%-------- Percent tip identity ####
ggslicepercid_bdtt <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>%
  ggplot() + geom_boxplot(aes(x=factor(slice), y=1-percid)) +
  # geom_hline(aes(yintercept=1-0.99), col="red") +
  xlab("Collapse depth into phylogenetic tree") +
  ylab("Percent difference between collapsed tips\n(log10 scale)") + scale_y_log10()  +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggslicepercid_bdtt
ggsave(filename = "07_FigureGeneration/percent_dissimilarity_across_slices_bdtt.png"
       ,ggslicepercid_bdtt)

#### %-------- Number of OTUs collapsed at each ####
gg_nnodes_collapsed <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>%
  ggplot() + geom_boxplot(aes(x=factor(slice), y=ncollapsed)) +
  scale_y_log10() + xlab("Collapse depth of phylogentic tree")+
  ylab("Number of tips collapsed into a node\n(log10 scale)")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  gg_nnodes_collapsed
ggsave(filename = "07_FigureGeneration/slice_depth_nNodes_collapsed.png"
       ,gg_nnodes_collapsed)

####%-------- Text summaries of collapsing ####
tab_ncollapsed_summary <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>% group_by(slice) %>%
  summarise(mean=mean(ncollapsed), geoMean=exp(mean(log(ncollapsed))))
write.table(tab_ncollapsed_summary, file="07_FigureGeneration/slices_ncollapsed_summary.txt", quote=FALSE, row.names = FALSE, sep="\t")


tab_slice_percid <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>% group_by(slice) %>%
  summarise(mean=mean(percid), geoMean=exp(mean(log(percid))))
write.table(tab_slice_percid, file="07_FigureGeneration/slices_PercID_summary.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### %-------- Taxonomic level collapse ####
gg_taxcollapse_bdtt <- allSlicesStats_bdtt %>%
  filter(ncollapsed>1) %>%
  mutate(taxlevel=ifelse(taxlevel==1,"Kingdom",ifelse(taxlevel==2, "Phylum"
                                                      , ifelse(taxlevel==3,"Class"
                                                               , ifelse(taxlevel==4,"Order"
                                                                        , ifelse(taxlevel==5,"Family"
                                                                                 ,ifelse(taxlevel==6,"Genus","Species"))))))) %>%
  mutate(taxlevel=factor(taxlevel, levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
  mutate(pres=1) %>%
  group_by(slice, taxlevel) %>%
  summarise(count=sum(pres)) %>% ungroup() %>%
  group_by(slice) %>% mutate(relabund=count/sum(count)) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x=factor(slice),y=relabund, fill=factor(taxlevel)), stat="identity") +
  xlab("Collapse depth into phylogenetic tree") +
  ylab("Proportion of nodes")+
  labs(fill="Highest resolution in taxonomy\nshared by all collapsed nodes") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "top"
        , legend.key.size = unit(0.5,"cm"), legend.text=element_text(size=7), legend.title = element_text(size=9))
gg_taxcollapse_bdtt
ggsave(filename = "07_FigureGeneration/taxonomy_collapsed_byslice_bdtt.png", width=6, height=4
       ,gg_taxcollapse_bdtt)


#### %-------- Combine percent and taxonomic collapse ####
gg_bdttcollapse <- plot_grid(ggslicepercid_bdtt+labs(title="a."),gg_taxcollapse_bdtt+labs(title="b."), ncol=2, align="h")
gg_bdttcollapse
ggsave(filename = "07_FigureGeneration/collapsing_bdtt_twopanel.png",width=10, height=5
       ,gg_bdttcollapse)

##### Metabolic function ########
#### %-------- FAPRO most common ####

fapro <- read.delim("04c_FAPROTAX/func_table.tsv") 
fapro <- fapro[, c("X.group", (meta[meta$sampleid %in% colnames(fapro),"sampleid"]))]
rnfapro <- fapro$X.group
fapro <- fapro[,-which(colnames(fapro)=="X.group")]
rownames(fapro)=rnfapro

fapro_withMeta <- fapro %>%
  as.data.frame() %>% rownames_to_column(var="MetabolicFunction") %>%
  pivot_longer(-MetabolicFunction, names_to="sampleid", values_to="count") %>%
  left_join(meta) 
gg_most_common_fapro <- fapro_withMeta %>%
  group_by(MetabolicFunction) %>% mutate(meanInFunc=sum(count)) %>% ungroup() %>%
  filter(meanInFunc>0) %>%
  arrange(-meanInFunc) %>% mutate(MetabolicFunction=factor(MetabolicFunction, levels=unique(MetabolicFunction))) %>%
  mutate(countLog=log10(count+1)) %>%
  ggplot() + geom_boxplot(aes(x=MetabolicFunction, y=countLog))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Counts observed per sample\n(log10)")
gg_most_common_fapro
ggsave(filename = "07_FigureGeneration/fapro_most_common.png", height=6, width=8
       , gg_most_common_fapro)

#### %-------- Metagenomic separate R2 ####

load("05a_metagenomic_inference/allMetaCor_bc.RData")
load("05a_metagenomic_inference/allMetaCor_pooled_bc.RData")
pooledOnly <-allMetaCor_pooled_bc %>%
  filter(metabType!="fapro") %>%
  mutate(predictor="pooled") %>%
  select(predictor, R2, pval, metabType) %>% distinct()

gg_meta_separateR2 <- allMetaCor_bc %>%
  filter(metabType!="fapro") %>%
  mutate(predictor=paste0(type,"_filt"))%>%
  full_join(pooledOnly) %>%
  left_join(mapPred) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval,n=5)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  mutate(metabType=ifelse(metabType=="ec","EC database","KO database")) %>%
  arrange(-R2) %>%
  mutate(Predictor=factor(Predictor, levels=unique(Predictor))) %>%
  ggplot(aes(x=Predictor, y=R2, fill=Predictor)) + 
  geom_bar(stat="identity", show.legend = FALSE) +
  facet_grid(metabType~.) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  # scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  scale_fill_manual(values= predColors) +
  ylab(expression(Variance~explained~by~predictor~(R^2))) +
  xlab("Predictor(s) used in regression") 
gg_meta_separateR2
ggsave(filename = "07_FigureGeneration/metab_separate_R2.png", height=8, width=4
       ,gg_meta_separateR2)


#### %-------- Text Metagenomic drop R2 calculations ####
maxR2 <- allMetaCor_pooled_bc %>% select(R2, metabType) %>%
  filter(metabType!="fapro") %>%
  rename(R2_total=R2) %>% distinct()
tab_metabR2drop <- allMetaCor_dropone %>%
  full_join(allMetaCor_pooled_bc) %>%
  filter(metabType!="fapro") %>%
  select(R2, type, metabType) %>% distinct() %>%
  left_join(maxR2) %>%
  mutate(R2_drop=R2_total-R2) %>%
  select(metabType, type, R2_drop, R2) %>%
  group_by(metabType) %>%
  arrange(-R2_drop, .by_group=TRUE)
tab_metabR2drop
write.table(tab_metabR2drop, file="07_FigureGeneration/table_metab_R2drop.txt", quote=FALSE, row.names = FALSE, sep="\t")

#### %-------- Metagenomic drop R2 ####
load("05a_metagenomic_inference/allMetaCor_dropone.RData")
# gg_metab_drop <- rbind(allMetaCor_dropone, allMetaCor_pooled_bc) %>%
#   filter(predictor!="Int", metabType!="fapro") %>%
#   select(R2, pval, type, metabType) %>%
#   rename(drop=type) %>%
#   mutate(drop=gsub("Dist","",gsub("drop_","no", drop))) %>%
#   mutate(drop=ifelse(drop=="pooled","fulldataset",drop))%>%
#   left_join(mapPred) %>%distinct() %>%
#   rowwise() %>% 
#   mutate(padj=p.adjust(pval,n=5)) %>% ungroup() %>%
#   mutate(Significant=padj<0.05) %>%
#   arrange(-R2) %>%
#   mutate(Model=factor(Model, levels=unique(Model))) %>%
#   mutate(metabType=ifelse(metabType=="ec","EC database","KO database")) %>%
#   ggplot(aes(x=Model, y=R2, fill=Model)) + 
#   geom_bar(stat="identity", show.legend = FALSE) +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
#   facet_grid(metabType~.) +
#   scale_fill_manual(values= dropColors) +
#   # scale_color_manual(values=c(`TRUE`="red",`FALSE`="black")) +
#   ylab(expression(Variance~explained~by~predictors~(R^2))) +
#   xlab("Model")
# gg_metab_drop
# ggsave(filename = "07_FigureGeneration/metab_dropR2.png", width=4, height=7
#        ,gg_metab_drop)


gg_metab_drop <- tab_metabR2drop%>%ungroup() %>%
  filter(type!="pooled") %>%
  mutate(type=gsub("Dist","",gsub("drop_","no", type))) %>%
  mutate(type=ifelse(type=="pooled","fulldataset",type))%>%
  rename(drop=type) %>%
  left_join(mapPred) %>%distinct() %>%
  mutate(R2_drop=-1*R2_drop) %>%
  arrange(R2_drop) %>%
  mutate(Model=factor(Model, levels=unique(Model))) %>%
  mutate(metabType=ifelse(metabType=="ec","EC database","KO database")) %>%
  ggplot(aes(x=Model, y=R2_drop, fill=Model)) + 
  geom_bar(stat="identity", show.legend = FALSE, alpha=0.5) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  facet_grid(metabType~.) +
  scale_fill_manual(values= dropColors) +
  # scale_color_manual(values=c(`TRUE`="red",`FALSE`="black")) +
  ylab(expression(Change~"in"~(R^2)~when~each~predictor~is~omitted)) +
  xlab("Model description")
gg_metab_drop
ggsave(filename = "07_FigureGeneration/metab_dropR2.png", width=4, height=7
       ,gg_metab_drop)
#### %-------- Combined metagenomic R2 ####
gg_aligned_meta <- plot_grid(gg_meta_separateR2 + labs(title="a."),gg_metab_drop+labs(title="b."), ncol=2, align="h")

ggsave(filename = "07_FigureGeneration/metab_two_panel_figure.png", height=7, width=8,
       gg_aligned_meta)

#### Prev/Excl ####

#### %-------- Prev/Excl OTUs ####
pe <- read.delim("06a_prevalence_and_exclusion/all_prev_and_excl_all_combos.txt")

gg_importantOTUs <- pe %>% 
  filter(PrevAndExclu_SpwithinLoc) %>%
  select(Host,OTU_fix,nSites) %>%
  mutate(Host=gsub("_"," ",Host)) %>%
  mutate(OTU_fix=ifelse(OTU_fix=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium.24 (Genus)", "Allo-Neo-Para-Rhizobium.24 (Genus)", OTU_fix)) %>%
  arrange(OTU_fix) %>%
  mutate(nSites=factor(nSites)) %>%
  ggplot() + geom_point(aes(y=OTU_fix, x=Host, cex=nSites)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  # scale_x_discrete(position = "top") +
  scale_size_manual(values=c(`1`=1,`2`=4))+
  labs(size="Number of sites\nOTU was both prevalent\nand exclusive on species") +
  xlab("Amphibian host species") +ylab("OTU")
gg_importantOTUs
ggsave(filename = "07_FigureGeneration/prevalent_and_exclusive_OTUs.png", height=8, width=8
       ,gg_importantOTUs)

##### BDTT vs CHEN #####

#### %-------- Chen coefficient ####

load("05b_BDTT_slicing/Bray_rare_chen/Correlations_bray_pooled_chen/allDat.RData")
allDatChen <- allDat
gg_chen_coef <- allDatChen %>%
  filter(type=="Full")%>%
  full_join(mapPred) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef,n=6)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  ggplot(aes(x=slice, y=abs(Coef), col=Predictor)) + 
  geom_line() +
  geom_point(aes(pch=Significant), cex=2) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),legend.key.size = unit(0.75,"cm"), legend.text = element_text(size=8))+
  scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  scale_color_manual(values= predColors) +
  ylab("Regression coefficient of predictor in pooled model") +
  xlab("Maximum depth of node for all tips\nfor all tips to be collapsed")
gg_chen_coef
ggsave(filename = "07_FigureGeneration/BDTT_chen_coef.png", height=5, width=7
       ,gg_chen_coef)

#### %-------- BDTT coefficient ####
gg_bdtt_coef <- allDat_bdtt_withpooled %>%
  filter(type=="Full")%>%
  full_join(mapPred) %>%
  rowwise() %>%
  mutate(padj=p.adjust(pval_coef,n=6)) %>% ungroup() %>%
  mutate(Significant=padj<0.05) %>%
  ggplot(aes(x=slice, y=abs(Coef), col=Predictor)) + 
  geom_line() +
  geom_point(aes(pch=Significant), cex=2) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),legend.key.size = unit(0.75,"cm"), legend.text = element_text(size=8))+
  scale_shape_manual(values=c(`TRUE`=19, `FALSE`=21)) +
  scale_color_manual(values= predColors) +
  ylab("Regression coefficient of predictor in pooled model") +
  xlab("Depth into phylogenetic tree\nat which tips are collapsed")
  gg_bdtt_coef
ggsave(filename = "07_FigureGeneration/BDTT_bdtt_coef.png", height=5, width=7
       ,gg_bdtt_coef)

#### %-------- Chen and BDTT combined ####
gg_aligned_bdttchen <- plot_grid(
  plot_grid(
    gg_chen_coef + labs(title="a.") + theme(legend.position = "none")
    , gg_bdtt_coef+labs(title="b.") + theme(legend.position = "none")
    , ncol = 2
    , align = "h")
  , plot_grid(
    get_legend(gg_chen_coef )
    , ncol =1)
  , rel_widths = c(7,2)
)

gg_aligned_bdttchen
ggsave(filename = "07_FigureGeneration/BDTT_two_panel_figure_chenvsbdtt.png", height=4, width=8,
       gg_aligned_bdttchen)
