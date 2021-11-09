#!bin/bash Rscript

#!bin/bash
library(picante) #loaded to use the example dataset
library(ape) # to get the Branch * sites matrices
library(abind)
library(caper) # to get the Branch * sites matrices
library(betapart) # to compute Beta-Diversity matrices
library(ecodist) # to compute correlations profiles (function 'MRM')
library(bigmemory)
library(castor)
library(phytools)
library(tidyverse)


dir.create("05b_BDTT_slicing_originalBDTT")

##### Load files ####

keepAbund <- read.delim("04e_partition_data_for_downstream/abundant_amphibians.txt")
meta <- read.delim("04b_data_filt_and_combine/downstream/meta_cut.txt") %>% 
  filter(Host %in% keepAbund$Host, Location %in% keepAbund$Location, CollectionDate %in% keepAbund$CollectionDate) %>%
  select(sampleid, Study, Host, Location, Month, HabitatClass, Latitude, Longitude, AmphTaxon, AmphFamily, AmphSubFamily)
# select(sampleid, Host, Study, Location, Month, CollectionDate, starts_with("biom_"), Latitude, Longitude, HabitatClass) 
otu <- read.delim("04b_data_filt_and_combine/downstream/otu_rare.txt", row.names=1) %>%
  select(one_of(meta$sampleid))
otu <- otu[rowSums(otu)>0,]

bacttree <- read.tree("04b_data_filt_and_combine/downstream/tree_filt_cal.nwk")
bacttree <- keep.tip(bacttree, rownames(otu))

print("Finished loading")
source("code/BDTT-master/BDTT_functions_MYC.R")
## RUNNING 
slices=round(seq(from=0, to=2, length.out=40),2) 
# 
# runBDTTOnly_withhighlow(slices=slices
#                         , intree=bacttree
#                         , sitesp = t(otu)
#                         , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen/"
#                         , sliceMethod = "bdtt"
#                         , doBeta = c("bc"))

# RUN ADDITIONAL SLICES
slices2=round(seq(from=0, to=0.3, length.out=31),2)
slices2 = slices2[!slices2%in%slices]
runBDTTOnly_withhighlow(slices=slices2
                        , intree=bacttree
                        , sitesp = t(otu)
                        , output = "05b_BDTT_slicing_originalBDTT/Bray_rare_chen_2/"
                        , sliceMethod = "bdtt"
                        , doBeta = c("bc")
                        , levels="Full"
                        )
