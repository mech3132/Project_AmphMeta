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

print("doing drop phylo")
predictorNamedList <- list(geoDist=geoDist, ecoDist=ecoDist, climDist=climDist, timeDist=timeDist, studyDist=studyDist)
slices_bdtt=round(seq(from=0, to=2, length.out=40),2)[1:10]
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_NoPhylo"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)
print("doing drop geo")
predictorNamedList <- list(phyloDist=phyloDist, ecoDist=ecoDist, climDist=climDist, timeDist=timeDist, studyDist=studyDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_NoGeo"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)

print("doing drop eco")
predictorNamedList <- list(geoDist=geoDist, phyloDist=phyloDist, climDist=climDist, timeDist=timeDist, studyDist=studyDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_NoEco"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)


print("doing drop clim")
predictorNamedList <- list(geoDist=geoDist, phyloDist=phyloDist, ecoDist=ecoDist, timeDist=timeDist, studyDist=studyDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_NoClim"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)

print("doing drop time")
predictorNamedList <- list(geoDist=geoDist, phyloDist=phyloDist, ecoDist=ecoDist, climDist=climDist, studyDist=studyDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_Notime"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)


print("doing drop study")
predictorNamedList <- list(geoDist=geoDist, phyloDist=phyloDist, ecoDist=ecoDist, climDist=climDist, timeDist=timeDist)
cors_bray_rare_chen <- runCorOnly_withPrevSlice_allpredictors(slices=slices_bdtt
                                                              , sites=meta$sampleid
                                                              , subsetName="bray_rare_NoStudy"
                                                              , output="05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
                                                              , metrics="bc"
                                                              , predictorNamedList=predictorNamedList
                                                              , levels=c("Full")
)
