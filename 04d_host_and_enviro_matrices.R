#!bin/bash
library(MASS)
library(ape)
library(vegan)
library(maps)
library(maptools)
library(ggfortify)
library(tidyverse)

dir.create("04d_host_and_enviro_matrices")
#### Preliminary analysis #####

# otu <- read.delim("04b_data_filt_and_combine/downstream/otu_rare.txt")
meta <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt")

amphtree <- read.tree("02c_generate_amphibian_phylogeny/amphibian_phylogeny.tre")
amphtree <- keep.tip(amphtree, unique(meta$Host))

#### Bioclimate ####
envFilt <- c("biom_AnnualMeanTemp", "biom_MeanDiurnalTempRange", "biom_TempSeasonality", "biom_TempAnnualRange", "biom_AnnualPrecip", "biom_PrecipSeasonality"  )
# meta_climate <- meta %>% select(Longitude, Latitude, starts_with("biom_")) %>% distinct() %>% drop_na() %>%
#   rename_at(vars(starts_with("biom_")), ~gsub("biom_","",.))
meta_climate <- meta %>% select(Longitude, Latitude, one_of(envFilt)) %>% distinct() %>% drop_na()
meta_climate_st <- apply(meta_climate[,-c(1,2)],2, function(x) (x-mean(x))/sd(x))
pca_env <- prcomp(meta_climate_st)
autoplot(pca_env, data=meta_climate
         , loadings=TRUE
         , loadings.label=TRUE)
ggsave("04d_host_and_enviro_matrices/pca_of_bioclim.png", height=3.5, width=4.5
       ,autoplot(pca_env, data=meta_climate
                 , loadings=TRUE
                 , loadings.label=TRUE)
       )
# So, let's make a combined distance metric of bioclimate data 
meta_bioclim_full <- apply(meta[,envFilt], 2, function(x) (x-mean(x, na.rm=TRUE)/sd(x, na.rm=TRUE)) )
rownames(meta_bioclim_full) <- meta$sampleid
climDist <- dist(meta_bioclim_full)


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
blueRamp <- colorRampPalette(c("hotpink","mediumblue"))
precipString <- latlong %>% left_join(meta %>% select(Latitude, Longitude, biom_AnnualPrecip)) %>% distinct() %>% pull(biom_AnnualPrecip) 
precipString_st <- (precipString-mean(precipString, na.rm=TRUE))/sd(precipString, na.rm=TRUE)
precipString_st <- precipString_st+abs(min(precipString_st, na.rm = TRUE))
# hist(precipString_st)
blueRampPrecip <- blueRamp(max(precipString_st, na.rm=TRUE)*100+1)
precipColString <- blueRampPrecip[round(precipString_st*100)+1]

# Get color ramp for richness
purpleRamp <- colorRampPalette(c("white","purple"))
aveRich <-meta %>% group_by(Latitude, Longitude) %>% summarise(meanRich=mean(log(observed_features), na.rm=TRUE))
richString <- latlong %>% left_join(aveRich) %>% distinct() %>% pull(meanRich)
purpleRampRich <- purpleRamp(max(richString, na.rm=TRUE))
richColString <- purpleRampRich[richString]# Get size ramp for number of samples
# Get number of samples
nSampDat <- meta %>% group_by(Latitude, Longitude) %>% summarise(nSamp=n()) %>% ungroup() %>%
  mutate(size=(nSamp+20)*10/max(nSamp+20)) 
nSampDat <- left_join(latlong, nSampDat)
# Sort so smaller points are on top
nSampDat$seasoncolor <- seasonColString
nSampDat$richcolor <- richColString
nSampDat$precipColor <- precipColString
nSampDat <- nSampDat %>% arrange(-size)
# Plot season
png("04d_host_and_enviro_matrices/map_samples_sizeSamp_colSeason.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$seasoncolor, pch=21, cex=nSampDat$size)
dev.off()
png("04d_host_and_enviro_matrices/map_samples_sizeSamp_colPrecip.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="white", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$precipColor, pch=21, cex=nSampDat$size)
dev.off()

png("04d_host_and_enviro_matrices/map_samples_colSeason.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$seasoncolor, pch=21, cex=1.5)
dev.off()
png("04d_host_and_enviro_matrices/map_samples_colPrecip.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$precipColor, pch=21, cex=1.5)
dev.off()

png("04d_host_and_enviro_matrices/map_samples_sizeSamp_colRich.png", height=600, width=900)
maps::map("world", fill=TRUE, col="lightgrey", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg=nSampDat$richcolor, pch=21, cex=nSampDat$size)
dev.off()

png("04d_host_and_enviro_matrices/map_samples_PLAIN.png", height=600, width=800)
maps::map("world", fill=TRUE, col="lightgrey", bg="white", ylim=c(-60, 90), mar=c(0,0,0,0))
points(nSampDat$Longitude,nSampDat$Latitude,bg="red", pch=21, cex=1.5)
dev.off()

# Make distance matrix
meta_geo <- meta %>% select(Latitude, Longitude)
rownames(meta_geo) <- meta$sampleid
geoDist <- dist(meta_geo)


### Amphibian Tree ####
strHab <- meta$HabitatClass[match(amphtree$tip.label, meta$Host)]
habcols <- data.frame(Habitat=c("Aquatic","Semi-aquatic","Arboreal","Terrestrial"), habCol=c("blue","purple","green","darkred"))
strHabCol <- left_join(data.frame(Habitat=strHab),habcols) %>% pull(habCol)
png("04d_host_and_enviro_matrices/amphibians_with_eco.png", height=1200, width = 800)
plot(amphtree, tip.color=strHabCol)
legend("topleft", c(habcols$Habitat), col= habcols$habCol, pch=19, cex=2)
dev.off()

png("04d_host_and_enviro_matrices/amphibians_plain.png", height=1200, width = 800)
plot(amphtree, tip.color="black")
# legend("topleft", c(habcols$Habitat), col= habcols$habCol, pch=19, cex=2)
dev.off()

# Make phyloDist
mattemp <- meta %>% select(sampleid, Host) %>% mutate(same=1) %>%
  pivot_wider(names_from=Host, values_from=same, values_fill = 0) %>%
  as.data.frame()
rownames(mattemp) <- mattemp[,'sampleid']
sitesp_mat <- as.matrix(mattemp[,-1])
amphdist <- as.matrix(cophenetic.phylo(amphtree))[colnames(sitesp_mat), colnames(sitesp_mat)]

# Multiply
phyloDist <- sitesp_mat%*%(amphdist %*%t(sitesp_mat))
# Check work
amphdist['Salamandra_salamandra','Ichthyosaura_alpestris']
meta %>% filter(sampleid %in% c("SRR7614958","MVB2015.2366")) %>%
  select(sampleid, Host)
phyloDist["SRR7614958","MVB2015.2366"]
phyloDist <- as.dist(phyloDist)

### Eco types ####
# define differences between ecotypes
listecotypes <- meta %>% select(HabitatClass) %>% drop_na() %>% pull() %>% unique()
ecoMat <- matrix(nrow=4, ncol=4, dimnames = list(listecotypes,listecotypes))
ecoMat[1,] <- c(0,1,1,2)
ecoMat[2,] <- c(1,0,1,2)
ecoMat[3,] <- c(1,1,0,1)
ecoMat[4,] <- c(2,2,1,0)
### Manually set ecoMat; Aquatic is more different than the other two
ecoTemp <- meta %>% select(sampleid, HabitatClass) %>% drop_na() %>%
  mutate(yes=1) %>%
  pivot_wider(names_from=HabitatClass, values_from=yes, values_fill = 0) %>%
  as.data.frame()
rownames(ecoTemp) <- ecoTemp$sampleid
siteEco_mat <- as.matrix(ecoTemp[,-1])
# Sort order
ecoMat <- ecoMat[colnames(siteEco_mat), colnames(siteEco_mat)]
# Multiply
ecoDist <- siteEco_mat %*%(ecoMat %*%t(siteEco_mat))

ecoDist <- as.dist(ecoDist)

#### timeDist ####
timeMeta <- meta %>% select(sampleid, Month) %>% drop_na() 
timeMetaMat <- as.matrix(timeMeta[,-1])
rownames(timeMetaMat) <- timeMeta$sampleid
timeDist <- as.matrix(dist(timeMetaMat))
timeDist[timeDist>6] <- 12-timeDist[timeDist>6]
# timeDist["J3.21","Ag01.2.NV"]
timeDist <- as.dist(timeDist)


#### funcDist ####

#### studyDist ####
studyonly <- meta %>% select(sampleid, Study) %>% mutate(yes=1) %>%
  pivot_wider(names_from=Study, values_from=yes, values_fill = 0) %>%
  as.data.frame()
rownames(studyonly) <- studyonly$sampleid
studyMat <- studyonly[,-1]
studyDist <- dist(studyMat)


#### Save ####
dir.create("04d_host_and_enviro_matrices/RData/")
save(ecoDist, file="04d_host_and_enviro_matrices/RData/ecoDist.RData")
save(phyloDist, file="04d_host_and_enviro_matrices/RData/phyloDist.RData")
save(climDist, file="04d_host_and_enviro_matrices/RData/climDist.RData")
save(geoDist, file="04d_host_and_enviro_matrices/RData/geoDist.RData")
save(timeDist, file="04d_host_and_enviro_matrices/RData/timeDist.RData")
save(studyDist, file="04d_host_and_enviro_matrices/RData/studyDist.RData")

# And also as tsv matrices
dir.create("04d_host_and_enviro_matrices/tsv")
write.table(as.matrix(ecoDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/ecoDist.txt", row.names=FALSE, quote = FALSE, sep="\t")
write.table(as.matrix(phyloDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/phyloDist.txt", row.names=FALSE, quote = FALSE, sep="\t")
write.table(as.matrix(climDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/climDist.txt", row.names=FALSE, quote = FALSE, sep="\t")
write.table(as.matrix(geoDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/geoDist.txt", row.names=FALSE, quote = FALSE, sep="\t")
write.table(as.matrix(timeDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/timeDist.txt", row.names=FALSE, quote = FALSE, sep="\t")
write.table(as.matrix(studyDist) %>% as.data.frame() %>%rownames_to_column(var="sampleid"), file="04d_host_and_enviro_matrices/tsv/studyDist.txt", row.names=FALSE, quote = FALSE, sep="\t")


