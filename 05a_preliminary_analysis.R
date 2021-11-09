#!bin/bash
library(MASS)
library(ape)
library(vegan)
library(maps)
library(maptools)
library(tidyverse)
library(ggfortify)

dir.create("05a_preliminary_analysis")
#### Preliminary analysis #####

# otu <- read.delim("04b_data_filt_and_combine/downstream/otu_rare.txt")
# keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt") 
amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
amphtree <- keep.tip(amphtree, unique(meta$Host))

#### Bioclimate ####
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

ggsave("05a_preliminary_analysis/pca_of_bioclim.png", height=6, width=8
       ,gg_climatepcoa
       )
# So, let's make a combined distance metric of bioclimate data 
load("04d_host_and_enviro_matrices/RData/climDist.RData")
# meta_bioclim_full <- apply(meta[,c("biom_TempSeasonality", "biom_AnnualPrecip")], 2, function(x) (x-mean(x, na.rm=TRUE)/sd(x, na.rm=TRUE)) )
# rownames(meta_bioclim_full) <- meta$sampleid
# climDist <- dist(meta_bioclim_full)


#### Geographic location #####
# Plot map
latlong <- meta[,c("Latitude","Longitude")] %>% distinct()
# Get color ramp for seasonality
redRamp <- colorRampPalette(c("yellow","red"))
seasonString <- latlong %>% left_join(meta %>% select(Latitude, Longitude, biom_TempSeasonality)) %>% distinct() %>% pull(biom_TempSeasonality)
seasonString_st <- (seasonString-mean(seasonString, na.rm=TRUE))/sd(seasonString, na.rm=TRUE)
seasonString_st <- seasonString_st+abs(min(seasonString_st, na.rm = TRUE))
# hist(seasonString_st)
redRampSeason <- redRamp(round(max(seasonString_st, na.rm=TRUE)*100)+1)
seasonColString <- redRampSeason[round(seasonString_st*100)+1]

# Get color ramp for annual precip
blueRamp <- colorRampPalette(c("lightgreen","darkblue"))
precipString <- latlong %>% left_join(meta %>% select(Latitude, Longitude, biom_AnnualPrecip)) %>% distinct() %>% pull(biom_AnnualPrecip) 
precipString_st <- (precipString-mean(precipString, na.rm=TRUE))/sd(precipString, na.rm=TRUE)
precipString_st <- precipString_st+abs(min(precipString_st, na.rm = TRUE))
# hist(precipString_st)
blueRampPrecip <- blueRamp(max(precipString_st, na.rm=TRUE)*100+1)
precipColString <- blueRampPrecip[round(precipString_st*100)+1]

# Get color ramp for richness
# purpleRamp <- colorRampPalette(c("white","purple"))
# aveRich <-meta %>% group_by(Latitude, Longitude) %>% summarise(meanRich=mean(log(observed_features), na.rm=TRUE))
# richString <- latlong %>% left_join(aveRich) %>% distinct() %>% pull(meanRich)
# purpleRampRich <- purpleRamp(max(richString, na.rm=TRUE))
# richColString <- purpleRampRich[richString]# Get size ramp for number of samples
# Get number of samples
nSampDat <- meta %>% group_by(Latitude, Longitude) %>% summarise(nSamp=n()) %>% ungroup() %>%
  mutate(size=(nSamp+20)*10/max(nSamp+20)) 
nSampDat <- left_join(latlong, nSampDat)
# Sort so smaller points are on top
nSampDat$seasoncolor <- seasonColString
# nSampDat$richcolor <- richColString
nSampDat$precipColor <- precipColString
nSampDat <- nSampDat %>% arrange(-size)
# Plot season
png("05a_preliminary_analysis/map_samples_sizeSamp_colSeason.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$seasoncolor, pch=21, cex=nSampDat$size)
dev.off()
png("05a_preliminary_analysis/map_samples_sizeSamp_colPrecip.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$precipColor, pch=21, cex=nSampDat$size)
dev.off()
# png("05a_preliminary_analysis/map_samples_sizeSamp_colRich.png", height=600, width=900)
# maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
# points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$richcolor, pch=21, cex=nSampDat$size)
# dev.off()

# Make distance matrix
load("04d_host_and_enviro_matrices/RData/geoDist.RData")
# meta_geo <- meta %>% select(Latitude, Longitude)
# rownames(meta_geo) <- meta$sampleid
# geoDist <- dist(meta_geo)

### Amphibian Tree ####
strHab <- meta$HabitatClass[match(amphtree$tip.label, meta$Host)]
habcols <- data.frame(Habitat=c("Aquatic","Semi-aquatic","Arboreal","Terrestrial"), habCol=c("blue","purple","green","darkred"))
strHabCol <- left_join(data.frame(Habitat=strHab),habcols) %>% pull(habCol)
png("05a_preliminary_analysis/amphibians_with_eco.png", height=1200, width = 800)
plot(amphtree, tip.color=strHabCol)
legend("topleft", c(habcols$Habitat), col= habcols$habCol, pch=19, cex=2)
dev.off()

# Make phyloDist
load("04d_host_and_enviro_matrices/RData/phyloDist.RData")
# mattemp <- meta %>% select(sampleid, Host) %>% mutate(same=1) %>%
#   pivot_wider(names_from=Host, values_from=same, values_fill = 0) %>%
#   as.data.frame()
# rownames(mattemp) <- mattemp[,'sampleid']
# sitesp_mat <- as.matrix(mattemp[,-1])
# amphdist <- as.matrix(cophenetic.phylo(amphtree))[colnames(sitesp_mat), colnames(sitesp_mat)]

# Multiply
# phyloDist <- sitesp_mat%*%(amphdist %*%t(sitesp_mat))
# Check work
# amphdist['Salamandra_salamandra','Ichthyosaura_alpestris']
# meta %>% filter(sampleid %in% c("SRR7614958","MVB2015.2366")) %>%
  # select(sampleid, Host)
# phyloDist["SRR7614958","MVB2015.2366"]

### Eco types ####
# define differences between ecotypes
load("04d_host_and_enviro_matrices/RData/ecoDist.RData")

# listecotypes <- meta %>% select(HabitatClass) %>% drop_na() %>% pull() %>% unique()
# ecoMat <- matrix(nrow=4, ncol=4, dimnames = list(listecotypes,listecotypes))
# ecoMat[1,] <- c(0,1,1,2)
# ecoMat[2,] <- c(1,0,1,2)
# ecoMat[3,] <- c(1,1,0,1)
# ecoMat[4,] <- c(2,2,1,0)
# ### Manually set ecoMat; Aquatic is more different than the other two
# ecoTemp <- meta %>% select(sampleid, HabitatClass) %>% drop_na() %>%
#   mutate(yes=1) %>%
#   pivot_wider(names_from=HabitatClass, values_from=yes, values_fill = 0) %>%
#   as.data.frame()
# rownames(ecoTemp) <- ecoTemp$sampleid
# siteEco_mat <- as.matrix(ecoTemp[,-1])
# # Sort order
# ecoMat <- ecoMat[colnames(siteEco_mat), colnames(siteEco_mat)]
# # Multiply
# ecoDist <- siteEco_mat %*%(ecoMat %*%t(siteEco_mat))
length(unique(meta$Host))
### Time of samples ####
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
  guides(colour="none")

ggsave(filename = "05a_preliminary_analysis/calendar_samples.png", height=4, width=6,
       gg_calendar)
load("04d_host_and_enviro_matrices/RData/timeDist.RData")


#### All predictors; test if they are correlated
allFactors <- c("phyloDist","ecoDist","geoDist","timeDist","climDist")

# ### NEED TO RUN STILL
# dir.create("05a_preliminary_analysis/predictor_corr")
# if ( ! file.exists("05a_preliminary_analysis/factorCorrelationsMantel.RData") ) {
#   factorCorrelationsMantel <- as.data.frame(matrix(ncol=4, nrow=0, dimnames=list(NULL, c("pval","R2","f1","f2"))))
#   for ( i in 1:(length(allFactors)-1) ) {
#     # i=1
#     for (j in (i+1):length(allFactors)) {
#       # j=2
#       f1 = as.matrix(get(allFactors[i]))
#       f2 = as.matrix(get(allFactors[j]))
#       sites = intersect(colnames(f1), colnames(f2))
#       # sites = intersect(sites, meta$sampleid[1:20]) # For troubleshooting function
#       # Check for NAs
#       f1na <- names(which(is.na(f1[,1])))
#       f2na <- names(which(is.na(f2[,1])))
#       sites = sites[!sites%in%c(f1na,f2na)]
#       f1 = as.dist(f1[sites,sites])
#       f2 = as.dist(f2[sites,sites])
#       # any(is.na(f2))
#       cor(f1,f2)
#       if ( ! file.exists(paste0("05a_preliminary_analysis/predictor_corr/mantel_",allFactors[i],"_",allFactors[j],".RData")) ) {
#         results = mantel(f1,f2, method="spearman")
#         assign(paste0("mantel_",allFactors[i],"_",allFactors[j]), results)
#         save(paste0("mantel_",allFactors[i],"_",allFactors[j]), file=paste0("05a_preliminary_analysis/predictor_corr/mantel_",allFactors[i],"_",allFactors[j],".RData") )
#       } else {
#         load(paste0("05a_preliminary_analysis/predictor_corr/mantel_",allFactors[i],"_",allFactors[j],".RData"))
#       results <- get(paste0("mantel_",allFactors[i],"_",allFactors[j]))
#         }
#       factorCorrelationsMantel <- rbind(factorCorrelationsMantel
#                                         ,data.frame(pval=results$signif, R2=results$statistic, f1=allFactors[i], f2=allFactors[j]))
#       
#     }
#   }
#   save(factorCorrelationsMantel, file = "05a_preliminary_analysis/factorCorrelationsMantel.RData")
#   
# } else {
#   load("05a_preliminary_analysis/factorCorrelationsMantel.RData")
# }

#### Richness and diversity patterns ####
# meta$observed_features
