#!bin/bash

###### Run gathering OTUs ###########

dir.create("05a_prev_and_exclus_calculations")

library(ape)
library(castor)
library(phytools)
library(tidyverse)

######## Loading #####
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt") 
otu <- read.delim("04b_data_filt_and_combine/downstream/otu_cut.txt", row.names=1) %>%
  select(one_of(meta$sampleid))
otu <- otu[rowSums(otu)>0,]
otu <- apply(otu,2,function(x) x/sum(x))
tree <- read.tree("04b_data_filt_and_combine/downstream/tree_filt_cal.nwk")
tree <- keep.tip(tree, rownames(otu))
taxa <- read.delim("04b_data_filt_and_combine/downstream/taxonomy.txt")

# taxa_g <- taxa %>% select(FeatureID, Genus) %>% mutate(Genus=gsub("[.].*$","",Genus))
source("code/prevalence_and_exclusion_functions.R")

dir.create(paste0("05a_prev_and_exclus_calculations/tables"))
dir.create(paste0("05a_prev_and_exclus_calculations/long_summaries"))
dir.create(paste0("05a_prev_and_exclus_calculations/RA_summaries"))
dir.create(paste0("05a_prev_and_exclus_calculations/summary_dat"))
dir.create(paste0("05a_prev_and_exclus_calculations/hostspecific_by_location"))


##### Get Bout and Host Mats ####
loc_with_multiple_dates <- meta %>% group_by(Location) %>%
  summarise(nDate = length(unique(CollectionDate))) %>% 
  filter(nDate>1) %>% pull(Location)
length(loc_with_multiple_dates)
loc_with_dates_to_keep <- meta %>% 
  group_by(Location, CollectionDate) %>%
  summarise(nSamps=n()) %>% group_by(Location) %>%
  mutate(max=nSamps==max(nSamps)) %>% ungroup() %>%
  filter(max) %>% select(Location, CollectionDate)
#Get list of combined Host/Loc that have more than 5 samples
sampsToKeep <- meta %>% inner_join(loc_with_dates_to_keep) %>%
  group_by(Host, Location, CollectionDate) %>%
  summarise(nsamps=n()) %>% filter(nsamps>=5) %>% ungroup() %>%
  select(Host, Location, CollectionDate)

meta_bout_wide <- meta %>% inner_join(sampsToKeep) %>%
  # unite(Location, CollectionDate, col=bout, remove=FALSE) %>%
  select(sampleid, Location) %>%
  mutate(bout=gsub("[:]|,|-","_",Location), pres=1) %>% 
  select(sampleid, bout, pres) %>%
  pivot_wider(names_from=bout, values_from=pres, values_fill=0) %>% as.data.frame()
meta_bout_wide_mat <- as.matrix(meta_bout_wide[,-1])
rownames(meta_bout_wide_mat) <- meta_bout_wide$sampleid

# Get wide dataset for species-- also filter by the datae/time restrictions
meta_Host_wide <- meta %>% inner_join(sampsToKeep) %>%
  select(sampleid, Host) %>%
  mutate(pres=1) %>% 
  pivot_wider(names_from=Host, values_from=pres, values_fill=0)
meta_Host_wide_mat <- as.matrix(meta_Host_wide[,-1])
rownames(meta_Host_wide_mat) <- meta_Host_wide$sampleid

## get HostLoc
meta_HostLoc_wide <- meta %>% inner_join(sampsToKeep) %>%
  mutate(bout=gsub("[:]|,|-","_",Location), pres=1) %>% 
  unite(Host, bout, col="HostLoc", sep="__") %>%
  select(sampleid, HostLoc) %>%
  mutate(pres=1) %>% 
  pivot_wider(names_from=HostLoc, values_from=pres, values_fill=0)
meta_HostLoc_wide_mat <- as.matrix(meta_HostLoc_wide[,-1])
rownames(meta_HostLoc_wide_mat) <- meta_HostLoc_wide$sampleid

## Finally, get multi host within locations
locationsMultipleHosts<- meta %>% inner_join(sampsToKeep) %>%
  group_by(Location) %>%
  summarise(nHosts=length(unique(Host))) %>%
  filter(nHosts>1) %>% select(Location) %>%
  inner_join(sampsToKeep)
meta_LocMultHosts <- meta %>% inner_join(locationsMultipleHosts) %>%
  select(sampleid, Host, Location)


slices=round(seq(0,2,length.out=40),2)
for ( SLICE in slices[c(1,2)]) {
  # SLICE=0
  print(SLICE)
  ##### Collapse tree ####
  if ( SLICE == 0 ) {
    otu_g <- otu
    rownames(otu_g) <- make.unique(taxa[match(rownames(otu_g), taxa$FeatureID),c("Genus")])
    otu_gPA <- as.matrix(otu_g)
    otu_gPA[otu_gPA>0] <- 1
  } else {
    tips_col_mat <- gatherOTUBy(tree, SLICE=SLICE)
    
    otu_g <- as.matrix(tips_col_mat)%*%as.matrix(otu)
    otu_gPA <- otu_g
    otu_gPA[otu_gPA>0] <- 1
  }
  #### Begin loops ######
  dir.create(paste0("05a_prev_and_exclus_calculations/Slice",SLICE))
  
  for (cutoff in c(0.8)) {
    # cutoff=0.9
    dothis=TRUE
    if ( dothis ) {#### Prev and Exclu calcuations ####
    prev_Host <- getPrevalence(otuPA_raw = otu_gPA, meta_mat_raw = meta_Host_wide_mat, cutoff=cutoff)
    prev_Loc <- getPrevalence(otuPA_raw = otu_gPA, meta_mat_raw = meta_bout_wide_mat, cutoff=cutoff)
    prev_HostLoc <- getPrevalence(otuPA_raw = otu_gPA, meta_mat_raw = meta_HostLoc_wide_mat, cutoff=cutoff)
    
    excl_Host <- getExclusive(otuPA_raw = otu_gPA, meta_mat_raw = meta_Host_wide_mat, cutoff=cutoff)
    excl_Loc <- getExclusive(otuPA_raw = otu_gPA, meta_mat_raw = meta_bout_wide_mat, cutoff=cutoff)
    excl_HostLoc <- getExclusive(otuPA_raw = otu_gPA, meta_mat_raw = meta_HostLoc_wide_mat, cutoff=cutoff)
    
    nOTU_Host <- getTotalNumberTaxa(otuPA_raw = otu_gPA, meta_mat_raw = meta_Host_wide_mat, Grouping="Host")
    nOTU_Loc <- getTotalNumberTaxa(otuPA_raw = otu_gPA, meta_mat_raw = meta_bout_wide_mat, Grouping="Location")
    nOTU_HostLoc <- getTotalNumberTaxa(otuPA_raw = otu_gPA, meta_mat_raw = meta_HostLoc_wide_mat, Grouping="HostLoc")

    ## Save with ASVs
    saveprevhost <- as.data.frame(t(prev_Host))%>%rownames_to_column(var="taxa")
    save(saveprevhost, file=paste0("05a_prev_and_exclus_calculations/tables/prev_Host_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
    saveprevloc <- as.data.frame(t(prev_Loc))%>%rownames_to_column(var="taxa")
    save(saveprevloc, file=paste0("05a_prev_and_exclus_calculations/tables/prev_Loc_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
    saveexclhost <- as.data.frame(excl_Host)%>%rownames_to_column(var="taxa")
    save(saveexclhost, file=paste0("05a_prev_and_exclus_calculations/tables/excl_Host_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
    saveexclloc <- as.data.frame(excl_Loc)%>%rownames_to_column(var="taxa")
    save(saveexclloc, file=paste0("05a_prev_and_exclus_calculations/tables/excl_Loc_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
    
      # Combine prevalence and plot
    prev_Host_long <- prev_Host %>% as.data.frame() %>% rownames_to_column(var="Host") %>% 
      pivot_longer(-Host, names_to="OTU", values_to="Prev") %>%
      # filter(Prev>0) %>% 
      left_join(data.frame(nSamps=colSums(meta_Host_wide_mat)) %>%rownames_to_column(var="Host"))%>% 
      left_join(nOTU_Host) %>%
      group_by(Host, nSamps, nOTUs) %>% summarise(nPrev=sum(Prev)) %>% mutate(Grouping="Host")%>% ungroup()
    prev_Loc_long <- prev_Loc %>% as.data.frame() %>% rownames_to_column(var="Location") %>%
      pivot_longer(-Location, names_to="OTU", values_to="Prev") %>%
      # filter(Prev>0) %>%
      left_join(data.frame(nSamps=colSums(meta_bout_wide_mat)) %>%rownames_to_column(var="Location"))%>% 
      left_join(nOTU_Loc) %>%
      group_by(Location, nSamps, nOTUs) %>% summarise(nPrev=sum(Prev)) %>% mutate(Grouping="Location")%>% ungroup()
    prev_HostLoc_long <- prev_HostLoc %>% as.data.frame() %>% rownames_to_column(var="HostLoc") %>%
      pivot_longer(-HostLoc, names_to="OTU", values_to="Prev") %>%
      # filter(Prev>0) %>%
      left_join(data.frame(nSamps=colSums(meta_HostLoc_wide_mat)) %>%rownames_to_column(var="HostLoc"))%>% 
      left_join(nOTU_HostLoc) %>%
      group_by(HostLoc, nSamps, nOTUs) %>% summarise(nPrev=sum(Prev)) %>% mutate(Grouping="HostLoc") %>%
      separate(HostLoc, into=c("Host","Location"), sep="__")%>% ungroup()
    prev_Host_Loc_HostLoc_long <- full_join(full_join(prev_Host_long, prev_Loc_long),prev_HostLoc_long)

    gg_prevcompareisons <- prev_Host_Loc_HostLoc_long %>%
      ggplot() + geom_point(aes(x=(nOTUs), y=(nPrev), col=Grouping)) +
      xlab("Number of samples in grouping") + ylab(paste0("Number of taxa prevalent\nat ",cutoff," threshold")) +
      labs(title=paste0("OTUs collapsed by slice depth ",SLICE))  +
      scale_x_log10() + scale_y_log10()
    gg_prevcompareisons
    ggsave(filename = paste0("05a_prev_and_exclus_calculations/Slice",SLICE,"/Prev_comparisons_slice",SLICE,"_cutoff",cutoff,".png"),
           gg_prevcompareisons)
    
    # Combine exclusivity and plot
    excl_Host_long <- t(excl_Host) %>% as.data.frame() %>% rownames_to_column(var="Host") %>% 
      pivot_longer(-Host, names_to="OTU", values_to="Excl") %>%
      # filter(Excl>0) %>% 
      left_join(data.frame(nSamps=colSums(meta_Host_wide_mat)) %>%rownames_to_column(var="Host"))%>% 
      left_join(nOTU_Host) %>%
      group_by(Host, nSamps, nOTUs) %>% summarise(nExcl=sum(Excl)) %>% mutate(Grouping="Host")
    excl_Loc_long <- t(excl_Loc) %>% as.data.frame() %>% rownames_to_column(var="Location") %>%
      pivot_longer(-Location, names_to="OTU", values_to="Excl") %>%
      # filter(Excl>0) %>%
      left_join(data.frame(nSamps=colSums(meta_bout_wide_mat)) %>%rownames_to_column(var="Location"))%>% 
      left_join(nOTU_Loc) %>%
      group_by(Location, nSamps, nOTUs) %>% summarise(nExcl=sum(Excl)) %>% mutate(Grouping="Location")
    excl_Host_Loc_long <- t(excl_HostLoc) %>% as.data.frame() %>% rownames_to_column(var="HostLoc") %>%
      pivot_longer(-HostLoc, names_to="OTU", values_to="Excl") %>%
      # filter(Excl>0) %>%
      left_join(data.frame(nSamps=colSums(meta_HostLoc_wide_mat)) %>%rownames_to_column(var="HostLoc"))%>% 
      left_join(nOTU_HostLoc) %>%
      group_by(HostLoc, nSamps, nOTUs) %>% summarise(nExcl=sum(Excl)) %>% mutate(Grouping="HostLoc") %>%
      separate(HostLoc, into=c("Host","Location"), sep="__")%>% ungroup()
    excl_Host_Loc_HostLoc_long<- full_join(full_join(excl_Host_long, excl_Loc_long),excl_Host_Loc_long)

    gg_exclcompareisons <- excl_Host_Loc_HostLoc_long %>%
      ggplot() + geom_point(aes(x=nSamps, y=nExcl, col=Grouping)) +
      xlab("Number of samples in grouping") + ylab(paste0("Number of taxa exclusive\nat ",cutoff," threshold")) +
      labs(title=paste0("OTUs collapsed by slice depth ",SLICE)) +
      scale_y_log10() + scale_x_log10()
    gg_exclcompareisons
    ggsave(filename = paste0("05a_prev_and_exclus_calculations/Slice",SLICE,"/Excl_comparisons_slice",SLICE,"_cutoff",cutoff,".png"),
           gg_exclcompareisons)
    
    
    ## Save
    all_Host_Loc_HostLoc_long <- full_join(excl_Host_Loc_HostLoc_long,prev_Host_Loc_HostLoc_long)
    save(all_Host_Loc_HostLoc_long, file=paste0("05a_prev_and_exclus_calculations/long_summaries/all_Host_Loc_HostLong_long_slice",SLICE,"_cutoff",cutoff,".RData"))
    
    #### Get abundance of taxa ####
    prev_Host_RA <- getAbundanceOfTaxa(mat_sitesp = prev_Host, meta_wide_mat = meta_Host_wide_mat, otu_g = otu_g
                                       , group = "Host", type="Prevalent", SLICE = SLICE, cutoff = cutoff)
    prev_Loc_RA <- getAbundanceOfTaxa(mat_sitesp = prev_Loc, meta_wide_mat = meta_bout_wide_mat, otu_g = otu_g
                                      , group = "Location", type="Prevalent", SLICE = SLICE, cutoff = cutoff)
    prev_HostLoc_RA <- getAbundanceOfTaxa(mat_sitesp = prev_HostLoc, meta_wide_mat = meta_HostLoc_wide_mat, otu_g = otu_g
                                      , group = "HostLoc", type="Prevalent", SLICE = SLICE, cutoff = cutoff)
    
    excl_Host_RA <- getAbundanceOfTaxa(mat_sitesp = t(excl_Host), meta_wide_mat = meta_Host_wide_mat, otu_g = otu_g
                                       , group = "Host", type="Exclusive", SLICE = SLICE, cutoff = cutoff)
    excl_Loc_RA <- getAbundanceOfTaxa(mat_sitesp = t(excl_Loc), meta_wide_mat = meta_bout_wide_mat, otu_g = otu_g
                                      , group = "Location", type="Exclusive", SLICE = SLICE, cutoff = cutoff)
    excl_HostLoc_RA <- getAbundanceOfTaxa(mat_sitesp = t(excl_HostLoc), meta_wide_mat = meta_HostLoc_wide_mat, otu_g = otu_g
                                          , group = "HostLoc", type="Exclusive", SLICE = SLICE, cutoff = cutoff)
    
    all_prev_excl_host_loc_RA <- rbind(prev_Host_RA, prev_Loc_RA, excl_Host_RA, excl_Loc_RA,prev_HostLoc_RA,excl_HostLoc_RA)
    save(all_prev_excl_host_loc_RA, file=paste0("05a_prev_and_exclus_calculations/RA_summaries/RA_all_Host_Loc_long_slice",SLICE,"_cutoff",cutoff,".RData"))
    
    
    gg_relabund_prevexcl <- all_prev_excl_host_loc_RA %>%
      group_by(subGroup, Group, Type) %>% mutate(meanRA=mean(RA), medRA=median(RA)) %>%
      filter(meanRA>0) %>%
      ggplot() + geom_histogram(aes(x=log10(meanRA), fill=Group), position="identity", alpha=0.5, bins=50) +
      facet_grid(Type~., scales="free") +
      xlab("Mean relative abundance\nof taxa by group (log10)")
    gg_relabund_prevexcl
    ggsave(filename = paste0("05a_prev_and_exclus_calculations/Slice",SLICE,"/RelativeAbund_comparisons_slice",SLICE,"_cutoff",cutoff,".png"),
           gg_relabund_prevexcl)
    
    #### Summarize ####
    prev_dat <- prev_Host_Loc_HostLoc_long %>% 
      mutate(percOTUlogPrev=log10(nPrev+1)/log10(nOTUs+1)) %>%
      mutate(incidOTUlogPrev=log10(nPrev+1)/log10(nSamps+1)) %>%
      # mutate(percOTUlogPrev=(nPrev)/(nOTUs)) %>%
      # mutate(incidOTUlogPrev=(nPrev)/(nSamps)) %>%
      select(Grouping, percOTUlogPrev,incidOTUlogPrev)
    
    excl_dat <- excl_Host_Loc_HostLoc_long %>%
      mutate(percOTUlogExcl=log(nExcl+1)/log(nOTUs+1)) %>%
      mutate(incidOTUlogExcl=log(nExcl+1)/log(nSamps+1)) %>%
      select(Grouping, percOTUlogExcl,incidOTUlogExcl )
    
    save(prev_dat, file=paste0("05a_prev_and_exclus_calculations/summary_dat/prev_dat_slice",SLICE,"_cutoff",cutoff,".RData"))
    save(excl_dat, file=paste0("05a_prev_and_exclus_calculations/summary_dat/excl_dat_slice",SLICE,"_cutoff",cutoff,".RData"))
    }
    ######## Host specificity BY LOCATION #########
    
    allDat_host_in_loc <- data.frame()
    for ( l in unique(meta_LocMultHosts$Location) ) {
      # get metadata filt
      metatemp <- meta_LocMultHosts[meta_LocMultHosts$Location==l,] 
      # Get meta wide mat
      metatemp_wide <- metatemp %>% select(sampleid,Host) %>%
        mutate(pres=1) %>% pivot_wider(names_from=Host, values_from=pres, values_fill=0) %>% as.data.frame()
      metatemp_wide_mat <- as.matrix(metatemp_wide[,-1])
      rownames(metatemp_wide_mat) <- metatemp_wide$sampleid
      # filter OTU table
      otuloc<- otu_g[,metatemp[,"sampleid"]]  
      otuloc <- otuloc[rowSums(otuloc)>0,]
      otulocPA <- otuloc
      otulocPA[otulocPA>0] <- 1
      prevtemp <- getPrevalence(otuPA_raw = otulocPA, meta_mat_raw = metatemp_wide_mat, cutoff=cutoff)
      excltemp <- getExclusive(otuPA_raw = otulocPA, meta_mat_raw = metatemp_wide_mat, cutoff=cutoff)
      prevandexcl <- t(prevtemp)*(excltemp)
      abundExcl <- getAbundanceOfTaxa(mat_sitesp = t(excltemp), otu_g = otuloc, meta_wide_mat = metatemp_wide_mat, SLICE=SLICE, group = l,type = "Exclusive", cutoff = cutoff)
      abundPrev <- getAbundanceOfTaxa(mat_sitesp = prevtemp, otu_g = otuloc, meta_wide_mat = metatemp_wide_mat, SLICE=SLICE, group = l,type = "Prevalent", cutoff = cutoff)
      allDat_host_in_loc <- rbind(allDat_host_in_loc, full_join(abundExcl, abundPrev))
    }
    allDat_host_in_loc_filt <- allDat_host_in_loc %>%
      group_by(subGroup, Type) %>%
      mutate(nLoc=length(unique(Group))) %>%
      filter(nLoc>1) %>% ungroup()
    
    summary_all_hostspecific_by_location <- allDat_host_in_loc %>%# filter(subGroup == unique(allDat_host_in_loc$subGroup)[1], Type=="Exclusive") %>%
      group_by(subGroup, Type, OTU, Group) %>% summarise(meanRA=mean(RA),maxRA=max(RA)) %>% ungroup() %>%
      group_by(subGroup, Type, OTU) %>% mutate(nAppearance=n())
    
    save(summary_all_hostspecific_by_location, file=paste0("05a_prev_and_exclus_calculations/hostspecific_by_location/summary_slice",SLICE,"_cutoff",cutoff,".RData"))
    
    
  }
  
}

########## Restrict by Location #############
# First, analyze BY LOCATION, each location has multiple species
# Then, look at exclusivity of ASVs between species; how many ASVs do they actually differ by?

# Verify that indicator ASVs on species are NOT the same across sites

# First, filter dataset to locations that have multiple species


