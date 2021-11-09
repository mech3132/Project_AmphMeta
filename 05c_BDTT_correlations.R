#!bin/bash
library(abind)
library(ecodist) # to compute correlations profiles (function 'MRM')
library(tidyverse)

source("code/BDTT-master/BDTT_functions_MYC.R")

### Load phylos ####

# keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
# meta2 <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt") %>% 
#   filter(Host %in% keepAbund$Host, Location %in% keepAbund$Location, CollectionDate %in% keepAbund$CollectionDate) %>%
#   select(sampleid, Study, Host, Location, Month, HabitatClass, Latitude, Longitude, AmphTaxon, AmphFamily, AmphSubFamily)
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt")
load("04d_host_and_enviro_matrices/RData/climDist.RData")
load("04d_host_and_enviro_matrices/RData/ecoDist.RData")
load("04d_host_and_enviro_matrices/RData/geoDist.RData")
load("04d_host_and_enviro_matrices/RData/phyloDist.RData")
load("04d_host_and_enviro_matrices/RData/timeDist.RData")
load("04d_host_and_enviro_matrices/RData/studyDist.RData")

sites <- meta[sample(nrow(meta),size = 20),"sampleid"]

slices_bdtt=round(seq(from=0, to=2, length.out=40),2)

#### This is obsolete, because I am doing each one separately below
# cors_bray_rare_chen <- runCorOnly_withPrevSlice(slices=slices_bdtt
#                                                 , sites=meta$sampleid
#                                                 , subsetName="bray_rare_chen"
#                                                 , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
#                                                 , metrics="bc"
#                                                 , phyloDist=phyloDist
#                                                 , geoDist=geoDist
#                                                 , ecoDist=ecoDist
#                                                 , climDist=climDist
#                                                 , timeDist=timeDist
#                                                 , studyDist=studyDist)
print("doing phylo")
cors_phylo <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_chen_phyloFull"
                                                              , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics = "bc"
                                                              , predictorNamedList = list(phyloDist=phyloDist)
                                                              , levels = "Full"
)

print("doing geo")

cors_geo <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                     , sites=meta$sampleid
                                                     , subsetName="bray_rare_chen_geoFull"
                                                     , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                     , metrics = "bc"
                                                     , predictorNamedList = list(geoDist=geoDist)
                                                     , levels = "Full"
)

print("doing clim")

cors_clim <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                     , sites=meta$sampleid
                                                     , subsetName="bray_rare_chen_climFull"
                                                     , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                     , metrics = "bc"
                                                     , predictorNamedList = list(climDist=climDist)
                                                     , levels = "Full"
)
print("doing time")

cors_time <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                     , sites=meta$sampleid
                                                     , subsetName="bray_rare_chen_timeFull"
                                                     , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                     , metrics = "bc"
                                                     , predictorNamedList = list(timeDist=timeDist)
                                                     , levels = "Full"
)

print("doing eco")

cors_eco <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                     , sites=meta$sampleid
                                                     , subsetName="bray_rare_chen_ecoFull"
                                                     , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                     , metrics = "bc"
                                                     , predictorNamedList = list(ecoDist=ecoDist)
                                                     , levels = "Full"
)

print("doing study")

cors_study <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                     , sites=meta$sampleid
                                                     , subsetName="bray_rare_chen_studyFull"
                                                     , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                     , metrics = "bc"
                                                     , predictorNamedList = list(studyDist=studyDist)
                                                     , levels = "Full"
)

#### POOLED ANALYSIS ######
print("doing pooled analysis")

predictorNamedList <- list(phyloDist=phyloDist, geoDist=geoDist, ecoDist=ecoDist, climDist=climDist, timeDist=timeDist, studyDist=studyDist)
slices_bdtt=round(seq(from=0, to=2, length.out=40),2) 
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                , sites=meta$sampleid
                                                , subsetName="bray_rare_ALLTOGETHER"
                                                , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                , metrics="bc"
                                                , predictorNamedList=predictorNamedList
                                                )


##### Non-bdtt version ####
print("doing non-bdtt version, only R2 output")

slices=round(seq(from=0, to=1.5, length.out=40),2)
cors_bray_rare_chen <- runCorOnly_withPrevSlice(slices=slices
                                                , sites=meta$sampleid
                                                , subsetName="bray_rare_chen"
                                                , output="05b_BDTT_slicing/Bray_rare_chen/"
                                                , metrics="bc"
                                                , phyloDist=phyloDist
                                                , geoDist=geoDist
                                                , ecoDist=ecoDist
                                                , climDist=climDist
                                                , timeDist=timeDist
                                                , studyDist=studyDist)

print("doing non-bdtt version, full MRM output")
predictorNamedList <- list(phyloDist=phyloDist, geoDist=geoDist, ecoDist=ecoDist, climDist=climDist, timeDist=timeDist, studyDist=studyDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_pooled_chen"
                                                              , output="05b_BDTT_slicing/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
)

