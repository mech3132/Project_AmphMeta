#!bin/bash
# library(MASS)
# library(vegan)

library(tidyverse)

dir.create("06a_prevalence_and_exclusion")
dir.create("06a_prevalence_and_exclusion/are_exclusive_or_prev_taxa_found_across_sites")
dir.create("06a_prevalence_and_exclusion/hist_exclusiveorprev")
dir.create("06a_prevalence_and_exclusion/tab_exclu_and_prev")
meta <- read.delim("04e_partition_data_for_downstream/final_meta.txt")
# Remove Hosts that are only from one location
hostToKeep <- meta %>% select(Host, Location) %>% distinct() %>% group_by(Host) %>%
  summarise(nsites=n()) %>% filter(nsites>1) %>% pull(Host)
# Remove Location that have only one host
locToKeep <- meta %>% select(Host, Location) %>% distinct() %>% group_by(Location) %>%
  summarise(nhost=n()) %>% mutate(Location=gsub(":|-|,","_",Location)) %>% 
  filter(nhost>1) %>% pull(Location)

load(paste0("05a_prev_and_exclus_calculations/"))
round(seq(0,2.0,length.out=40),2)

slices<- c(0,0.05,0.26,0.36)
cutoffs <- c(0.8)

SLICE=slices[2]
cutoff=cutoffs[1]

################Looking at prevalence and exclusivity################
load(paste0("05a_prev_and_exclus_calculations/long_summaries/all_Host_Loc_HostLong_long_slice",SLICE,"_cutoff",cutoff,".RData"))
all_Host_Loc_HostLoc_long %>% ungroup() %>%
  filter(Grouping=="Host") %>% 
  mutate(percPrev=nPrev/nOTUs, prevPerSamp=nPrev/nSamps) %>%
  select(Host, nPrev, percPrev, prevPerSamp) %>%
  filter(Host %in% hostToKeep) %>%
  ggplot() + geom_histogram(aes(x=log10(percPrev)))

host_prevalence_summary <- all_Host_Loc_HostLoc_long %>% ungroup() %>%
  filter(Grouping=="Host") %>% 
  filter(Host %in% hostToKeep) %>%
  mutate(percPrev=nPrev/nOTUs, percExcl=nExcl/nOTUs) %>%
  select(Host,nSamps, nOTUs, nPrev, nExcl, percPrev, percExcl) 

########Prevalence
host_prevalence_summary_mat <- host_prevalence_summary%>% distinct() %>% 
  mutate(HostCount=1) %>%select(-Host) %>% as.matrix()
sink(file="06a_prevalence_and_exclusion/table_of_host_prevalence.txt")
apply(host_prevalence_summary_mat, 2, function(x) c(TOTAL=sum(x), min=min(x), max=max(x)
                                                    , geomean=exp(mean(log(x+1)))-1
                                                    , geoLower=exp(mean(log(x+1))-2*sd(log(x+1)))-1
                                                    , geoUpper=exp(mean(log(x+1))+2*sd(log(x+1)))-1
                                                    , mean=mean(x)))
sink()


location_prevalence_summary <- all_Host_Loc_HostLoc_long %>% ungroup() %>%
  filter(Grouping=="Location") %>% 
  filter(Location %in% locToKeep) %>%
  mutate(percPrev=nPrev/nOTUs, percExcl=nExcl/nOTUs) %>%
  select(Location,nSamps, nOTUs, nPrev, nExcl, percPrev, percExcl)
location_prevalence_summary_mat <- location_prevalence_summary%>% distinct() %>% 
  mutate(LocCount=1) %>%select(-Location) %>% as.matrix()
sink(file="06a_prevalence_and_exclusion/table_of_location_prevalence.txt")
apply(location_prevalence_summary_mat, 2, function(x) c(TOTAL=sum(x), min=min(x), max=max(x)
                                                        , geomean=exp(mean(log(x+1)))-1
                                                        , geoLower=exp(mean(log(x+1))-2*sd(log(x+1)))-1
                                                        , geoUpper=exp(mean(log(x+1))+2*sd(log(x+1)))-1
                                                        , mean=mean(x)))
sink()

## See if number of samples within group differs
test_nsamp_diff_hostloc <- all_Host_Loc_HostLoc_long %>% ungroup() %>%filter(Grouping %in% c("Host","Location")) %>%
  filter( (Grouping=="Host" & Host %in% hostToKeep) | (Grouping=="Location" & Location %in% locToKeep ) ) %>%
  select(Grouping, nSamps)
sink(file="06a_prevalence_and_exclusion/anova_nSamps_HostvsLoc.txt")
anova(lm(nSamps ~ Grouping, data=test_nsamp_diff_hostloc))
sink()

####  Prevalence ####

test_nPrev_diff_hostloc <- all_Host_Loc_HostLoc_long %>% ungroup() %>%filter(Grouping %in% c("Host","Location")) %>%
  filter( (Grouping=="Host" & Host %in% hostToKeep) | (Grouping=="Location" & Location %in% locToKeep ) ) %>%
  mutate(percPrev=nPrev/nOTUs) %>% mutate(logpercPrev=log(percPrev+1)) %>% 
  select(Grouping, logpercPrev)
sink(file="06a_prevalence_and_exclusion/anova_logpercPrev_HostvsLoc.txt")
anova(lm(logpercPrev ~ Grouping, data=test_nPrev_diff_hostloc))
sink()
#### Exclusive ####

test_nExcl_diff_hostloc <- all_Host_Loc_HostLoc_long %>% ungroup() %>%filter(Grouping %in% c("Host","Location")) %>%
  filter( (Grouping=="Host" & Host %in% hostToKeep) | (Grouping=="Location" & Location %in% locToKeep ) ) %>%
  mutate(percExcl=nExcl/nOTUs) %>% mutate(logpercExcl=log(percExcl+1)) %>% 
  select(Grouping, logpercExcl)
sink(file="06a_prevalence_and_exclusion/anova_logpercExcl_HostvsLoc.txt")
anova(lm(logpercExcl ~ Grouping, data=test_nExcl_diff_hostloc))
sink()
########### Prevalent AND exclusive acorss host and location? ##########
# Exclu and prev by host
load(paste0("05a_prev_and_exclus_calculations/tables/excl_Host_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
load(paste0("05a_prev_and_exclus_calculations/tables/prev_Host_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
host_excl_long <- saveexclhost %>%
  pivot_longer(-taxa, names_to="Host", values_to="Exclusive")
host_prev_long <- saveprevhost %>%
  pivot_longer(-taxa, names_to="Host", values_to="Prevalent")

# Combined
bothExcluAndPrev_HOST <- full_join(host_excl_long, host_prev_long) %>%
  filter(Host %in% hostToKeep) %>%
  mutate(Both=Exclusive*Prevalent) %>%
  filter(Both>0)
bothExcluAndPrev_HOST %>% arrange(Host)
### NONE ###

# Excl and prev by location 
load(paste0("05a_prev_and_exclus_calculations/tables/excl_Loc_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
load(paste0("05a_prev_and_exclus_calculations/tables/prev_Loc_mat_slice",SLICE,"_cutoff",cutoff,".RData"))
loc_excl_long <- saveexclloc %>%
  pivot_longer(-taxa, names_to="Location", values_to="Exclusive")
loc_prev_long <- saveprevloc %>%
  pivot_longer(-taxa, names_to="Location", values_to="Prevalent")

bothExcluAndPrev_LOCATION <- full_join(loc_excl_long, loc_prev_long) %>%
  filter(Location %in% locToKeep) %>%
  mutate(Both=Exclusive*Prevalent) %>%
  filter(Both>0)
bothExcluAndPrev_LOCATION %>% arrange(Location)

################ Looking at relative abundances################
load(paste0("05a_prev_and_exclus_calculations/RA_summaries/RA_all_Host_Loc_long_slice",SLICE,"_cutoff",cutoff,".RData"))
# all_prev_excl_host_loc_RA

host_RA_string <- all_prev_excl_host_loc_RA %>% filter(Group=="Host", Type=="Prevalent") %>% group_by(OTU) %>%
  summarise(meanRA=mean(log(RA+1))) %>% pull(meanRA)
loc_RA_string <- all_prev_excl_host_loc_RA %>% filter(Group=="Location", Type=="Prevalent") %>%  group_by(OTU) %>%
  summarise(meanRA=mean(log(RA+1))) %>% pull(meanRA)
ttest_RA_hostvsloc <- t.test(x=host_RA_string, y=loc_RA_string)
ttest_RA_hostvsloc

host_RAvar_string <- all_prev_excl_host_loc_RA %>% filter(Group=="Host", Type=="Prevalent") %>% group_by(OTU) %>%
  summarise(varRA=var(log(RA+1))) %>% pull(varRA)
loc_RAvar_string <- all_prev_excl_host_loc_RA %>% filter(Group=="Location", Type=="Prevalent") %>%  group_by(OTU) %>%
  summarise(varRA=var(log(RA+1))) %>% pull(varRA)
ttest_RA_hostvsloc <- t.test(x=host_RAvar_string, y=loc_RAvar_string)
ttest_RA_hostvsloc


########## Look at site-specific host patterns ########

load(paste0("05a_prev_and_exclus_calculations/hostspecific_by_location/summary_slice",SLICE,"_cutoff",cutoff,".RData"))

## Things that are prevalent tend to be abundant, while things that are exclusive tend to be low abundance.
gg_hist_exclu_or_prev <- summary_all_hostspecific_by_location %>% ungroup() %>%
  ggplot() + geom_histogram(aes(x=log10(meanRA), fill=Type), position="identity",  alpha=0.5) 
gg_hist_exclu_or_prev
ggsave(filename = paste0("06a_prevalence_and_exclusion/hist_exclusiveorprev/hist_slice",SLICE,"_cutoff",cutoff,".png")
       ,gg_hist_exclu_or_prev)

gg_hist_exclu_or_prev_filt <- summary_all_hostspecific_by_location %>% ungroup() %>%
  filter(maxRA>0.01) %>%
  ggplot() + geom_histogram(aes(x=log10(meanRA), fill=Type), position="identity",  alpha=0.5) 
gg_hist_exclu_or_prev_filt
ggsave(filename = paste0("06a_prevalence_and_exclusion/hist_exclusiveorprev/hist_FILTERED_slice",SLICE,"_cutoff",cutoff,".png")
       ,gg_hist_exclu_or_prev_filt)

gg_hist_exclu_or_prev_filt_maxRA <- summary_all_hostspecific_by_location %>% ungroup() %>%
  filter(maxRA>0.01) %>%
  ggplot() + geom_histogram(aes(x=log10(maxRA), fill=Type), position="identity",  alpha=0.5) 
gg_hist_exclu_or_prev_filt_maxRA
ggsave(filename = paste0("06a_prevalence_and_exclusion/hist_exclusiveorprev/hist_FILTERED_maxRA_slice",SLICE,"_cutoff",cutoff,".png")
       ,gg_hist_exclu_or_prev_filt_maxRA)

# List of OTUs that are prevalent within site
prevList <- summary_all_hostspecific_by_location %>% ungroup() %>%
  rename(nAppearance_prev=nAppearance) %>%
  filter(Type=="Prevalent") %>% select(subGroup,OTU,Group,nAppearance_prev) %>%
  mutate(Prevalent=TRUE)
exclList <- summary_all_hostspecific_by_location %>% ungroup() %>%
  filter(subGroup %in% hostToKeep) %>% # Filter to only species that are found with other species
  rename(nAppearance_excl=nAppearance) %>%
  filter(Type=="Exclusive") %>% select(subGroup,OTU,Group,nAppearance_excl) %>%
  mutate(Exclusive=TRUE)

full_summary_prev_and_excl_spwithinsite <- full_join(prevList, exclList) %>%
  mutate(nAppearance_excl=ifelse(is.na(nAppearance_excl), 0, nAppearance_excl)
         , nAppearance_prev=ifelse(is.na(nAppearance_prev),0,nAppearance_prev)
         ,Prevalent=ifelse(is.na(Prevalent),FALSE,Prevalent)
         , Exclusive=ifelse(is.na(Exclusive), FALSE, Exclusive)) %>%
  group_by(subGroup) %>% mutate(nSites=length(unique(Group))) %>% ungroup() 

# ASVs in each combination that are prevalent AND exclusive
toPrint <- full_summary_prev_and_excl_spwithinsite %>%
  filter(Prevalent & Exclusive) %>%
  select(subGroup, OTU, Group, Prevalent, Exclusive) %>%
  group_by(subGroup, OTU, Prevalent, Exclusive) %>% summarise(nSites=length(unique(Group))) %>%
  arrange(OTU)
toPrint
write.table(toPrint, file=paste0("06a_prevalence_and_exclusion/tab_exclu_and_prev/tab_slice",SLICE,"_cutoff",cutoff,".txt"),
            quote = FALSE, row.names = FALSE, sep="\t")


bothExcluAndPrev_HOST
bothExcluAndPrev_LOCATION

tocombine <- toPrint%>%ungroup() %>% rename(Host=subGroup) %>%
  select(Host, OTU, nSites) %>% mutate(PrevAndExclu_SpwithinLoc=TRUE)
tocombine_host <- bothExcluAndPrev_HOST %>% rename(OTU=taxa) %>% 
  select(Host, OTU) %>% mutate(PrevAndExclu_Host=TRUE)
tocombine_loc <- bothExcluAndPrev_LOCATION%>% rename(OTU=taxa) %>% 
  select(Location, OTU) %>% mutate(PrevAndExclu_Loc=TRUE)

allprevandeclu_alltrials <- full_join(tocombine, tocombine_host) %>% full_join(tocombine_loc) %>%
  select(Location, Host, OTU, everything())

listOTUsToFix <- allprevandeclu_alltrials$OTU
fixedOTUs <- data.frame(OTU=listOTUsToFix, OTU_fix=NA)
for ( o in listOTUsToFix ) {
  lvl <- gsub("__.*$","",o)
  suff <- ifelse(lvl=="g", "(Genus)", ifelse(lvl=="f","(Family)", ifelse(lvl=="o","(Order)", ifelse(lvl=="c","(Class)",NA))))
  pre <- gsub("^.*__","",o)
  fixedOTUs[which(fixedOTUs$OTU==o),"OTU_fix"] <- paste(pre,suff,sep=" ")
}
allprevandeclu_alltrials <- allprevandeclu_alltrials %>% left_join(fixedOTUs)

write.table(allprevandeclu_alltrials, file="06a_prevalence_and_exclusion/all_prev_and_excl_all_combos.txt", quote=FALSE, row.names = FALSE, sep="\t")

gg_importantOTUs <- allprevandeclu_alltrials %>% 
  filter(PrevAndExclu_SpwithinLoc) %>%
  select(Host,OTU_fix,nSites) %>%
  arrange(OTU_fix) %>%
  ggplot() + geom_point(aes(x=OTU_fix, y=Host, cex=factor(nSites))) +
  theme(axis.text.x = element_text(angle=45, hjust=0)) +
  scale_x_discrete(position = "top") 
gg_importantOTUs
ggsave(filename = "06a_prevalence_and_exclusion/import_otus_prevandexcl_spwithinloc.png", height=6, width=12
       ,gg_importantOTUs)
#summary
sink(file = "06a_prevalence_and_exclusion/summary_prev_and_excl_spwithloc_counts.txt")
print("number of host species")
allprevandeclu_alltrials %>% 
  filter(PrevAndExclu_SpwithinLoc) %>%
  select(Host,OTU_fix,nSites) %>%
  pull(Host) %>% unique() %>% length()
print("number of OTUs")
allprevandeclu_alltrials %>% 
  filter(PrevAndExclu_SpwithinLoc) %>%
  select(Host,OTU_fix,nSites) %>%
  pull(OTU_fix) %>% unique() %>% length()
sink()

