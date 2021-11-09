#!bin/bash
library(stringr)
library(raster)
library(tidyverse)
# Get ecology for amphibians

### 
# jord <- read.delim("00_manual_metadata/Mapping_global_July_21.txt")
meta <- read.delim("01_parsed_meta/full_meta.txt", sep="\t")
allHabClass <- read.delim("00_manual_metadata/allHabitatClasses_cleaned.tsv")
dir.create("02b_bioclim_ecology")


# This is copied from 01_parsed-meta so should be the same
# allHabClass <- allHabClass %>% select(FrogSpecies, HabitatClass) %>% distinct() %>% 
#   filter(!is.na(FrogSpecies)) %>%
#   rowwise() %>% mutate(unknownsp = length(grep("(sp|spp)[.]$",FrogSpecies))>0) %>% ungroup() %>%
#   rowwise() %>% mutate(cfsp = length(grep("cf.",FrogSpecies))>0) %>% ungroup() %>%
#   filter(!unknownsp, !cfsp , !FrogSpecies %in% c("Pseudoeurycea sp. n. Mozotal 2"
#                                           ,"Bolitoglossa franklini lincolni"
#                                           ,"Craugastor"
#                                           ,"Pseudoeurycea sp. n. Mozotal 1"
#                                           ,"Scinax aromothylla"
#                                           ,"Scinax gr. ruber")) %>%
#   mutate(FrogSpecies = gsub(" x ", " ", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Agalychnis granulosa", "Hylomantis granulosa", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Aplastodiscus eugenoi", "Aplastodiscus eugenioi", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Aplastodiscus pervidis", "Aplastodiscus perviridis", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Bolitoglossa franklini nigroflavescens", "Bolitoglossa franklini", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Lithobates actesbeianus", "Rana catesbeiana", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Lithobates catesbeianus", "Rana catesbeiana", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Lithobates catesbeiana", "Rana catesbeiana", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Lithobates clamitans", "Rana clamitans", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Lithobates sylvaticus", "Rana sylvatica", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Pithecopus rhodei", "Phyllomedusa rohdei", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Scinax flavogutattus", "Scinax flavoguttatus", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Scinax flavogutattus", "Scinax flavoguttatus", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Boophis goudoti", "Boophis goudotii", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Anodonthyla boulengeri", "Anodonthyla boulengerii", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Cophyla grandis", "Platypelis grandis", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Boophis boehemi", "Boophis boehmei", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("R. cat", "Rana catesbeiana", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Mantidactylus cowanii small", "Mantidactylus cowanii", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Hylomantis lemur", "Agalychnis lemur", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("Eleutherodactylus caspari", "Eleutherodactylus casparii", FrogSpecies)) %>%
#   mutate(FrogSpecies = gsub("iii", "ii", FrogSpecies)) %>% # IN case there was added extra
#   rename(Host=FrogSpecies) %>%
#   filter(Host %in% meta$Host) %>%
#   select(Host, HabitatClass) %>% distinct() %>%
#   mutate(HabitatClass=ifelse(HabitatClass=="",NA,HabitatClass)) %>%
#   mutate(HabitatClass = gsub("arboreal","Arboreal", HabitatClass))

# ### FIND DUPLICATE ECOS, APPEND WHERE NEEDED
# allHabClass_adj <- allHabClass %>% filter(!is.na(HabitatClass)) %>% distinct() %>%mutate(listed=TRUE) %>%
#   pivot_wider(names_from=HabitatClass, values_from=listed, values_fill = FALSE) %>%
#   rowwise() %>%
#   mutate(multiple=sum(across(-Host))) %>% ungroup() %>%
#   mutate(HabitatClass=ifelse(Terrestrial, "Terrestrial", ifelse(Arboreal, "Arboreal", ifelse(Aquatic, "Aquatic", ifelse(`Semi-aquatic`, "Semi-aquatic", NA)))))

# allHabClassFinal <- allHabClass_adj %>% select(Host, HabitatClass)
#TEST
# meta_filt <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt")

allHabClass_adj <- allHabClass %>% mutate(Host=gsub("_"," ",Host))
allHabClassFinal <- meta %>% left_join(allHabClass_adj) %>% select(Host, HabitatClass) %>% distinct() %>% arrange(HabitatClass)

# Write table for val
write.table(allHabClassFinal, file="02b_bioclim_ecology/allHabitateClasses.txt", row.names = FALSE, quote = FALSE, sep="\t")

###### BIOM CLIM #######
meta <- meta %>%
  mutate(Longitude=as.numeric(str_trim(Longitude)), Latitude=as.numeric(str_trim(Latitude)))

#### Bioclim data and amphibian data ####
BioClimData <- getData("worldclim", var="bio", res=2.5)
# test2 <- getData("worldclim", var="bio", res=10)

allVar <- BioClimData[[c(seq(1,19))]]
names(allVar) <- c("AnnualMeanTemp","MeanDiurnalTempRange","DiurnalTempRange/TempAnnualRange"
                   ,"TempSeasonality","MaxTempWarmestMonth","MinTempCodestMonth"
                   , "TempAnnualRange","MeanTempWettestQuarter","MeanTempDriestQuarter"
                   , "MeanTempWarmestQuarter","MeanTempColdestQuarter", "AnnualPrecip"
                   ,"PrecipWettestMonth", "PrecipDriestMonth","PrecipSeasonality"
                   , "PrecipWettestQuarter","PrecipDriestQuarter","PrecipWarmestQuarter","PrecipColdestQuarter")
names(allVar) <- paste0("biom_",names(allVar))
coords <- data.frame(y=meta$Latitude, x=meta$Longitude) %>% distinct() %>% drop_na() %>%
  select(x,y)
# coords <- data.frame(x=meta$Longitude,y=meta$Latitude)
points <- SpatialPoints(coords, proj4string = allVar@crs)
climVals <- raster::extract(allVar, points)
meta_withbio <- meta %>% left_join(data.frame(climVals, Longitude=coords[,"x"], Latitude=coords[,"y"]))


### Write out
meta_all <- meta_withbio %>% full_join(allHabClassFinal) 

write.table(meta_all, file="02b_bioclim_ecology/meta_complete.txt", row.names = FALSE, quote = FALSE, sep="\t")
