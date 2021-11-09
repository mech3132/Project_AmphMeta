#!bin/bash
# Temp
library(tidyverse)
library(lubridate)
library(stringr)

dir.create("01_parsed_meta")

amph <- read.delim("00_all_sample_meta/all_sample_meta.txt", na.strings = c("NA","","na","unknown","<NA>"))
# Get metadata of pooled samples
walke2015 <- read.delim("00_manual_metadata/walke2015_10272_20180418-081219.txt", na.strings = c("NA","na","","<NA>"))
walke_study_only <- amph %>% filter(`ENA.STUDY`=="SRP062395") 
walke_study_full <- walke_study_only[,which(apply(walke_study_only, 2, function(x) any(!is.na(x))))] %>%
  select(-c(sample_type)) %>%
  full_join(walke2015[,which(apply(walke2015, 2, function(x) any(!is.na(x))))] %>% mutate(`ENA.STUDY`="SRP062395")) %>%
  select(-altitude) %>% mutate(host = host_common_name)

longo2017 <- read.delim("00_manual_metadata/Mapping_file_Longomoleco.txt", , na.strings = c("NA","na","","<NA>")) 
longo_study_only <- amph %>% filter(`ENA.STUDY`=="SRP104695") 
longo_study_full <- longo_study_only[,which(apply(longo_study_only, 2, function(x) any(!is.na(x))))] %>%
  full_join(longo2017[,which(apply(longo2017, 2, function(x) any(!is.na(x))))] %>% mutate(`ENA.STUDY`="SRP104695")) %>%
  mutate(host = Species)

ruthsatz2020 <- read.delim("00_manual_metadata/ruthsatz.txt", , na.strings = c("NA","na","","<NA>")) 
ruthsatz2020 <- ruthsatz2020 %>% rename(SUBMITTER_ID.label.Sample.name=X.SampleID)
ruthsatz_study_only <- amph %>% filter(`ENA.STUDY`=="SRP285287") 
ruthsatz_study_full <- ruthsatz_study_only[,which(apply(ruthsatz_study_only, 2, function(x) any(!is.na(x))))] %>%
  full_join(ruthsatz2020[,which(apply(ruthsatz2020, 2, function(x) any(!is.na(x))))] %>% mutate(`ENA.STUDY`="SRP285287")) %>%
  mutate(host=Species)

kuenExtra <- read.delim("00_manual_metadata/Mapping_global_July_21.txt")
# kuenExtra2 <- read.delim("00_manual_metadata/Kuen_updated_mapping.txt")
# 
# tempsummary <- read.delim("00_nonEBI_data/Kuen/summary_og_table.txt", sep = ":")
# actualSamples <- tempsummary$Num.samples[13:length(tempsummary$Num.samples)]
# 
# sum(kuenExtra$X.SampleID %in% actualSamples)
# sum(kuenExtra2$X.SampleID %in% actualSamples)
# sum(unique(c(kuenExtra$X.SampleID, kuenExtra2$X.SampleID)) %in% actualSamples)
# 
# kuenExtraExclusive <- kuenExtra$X.SampleID[which(!kuenExtra$X.SampleID %in% kuenExtra2$X.SampleID)]
# sum(kuenExtraExclusive %in% actualSamples)
# 
# kuenExtra2$X.SampleID[which(!kuenExtra2$X.SampleID %in% kuenExtra$X.SampleID)]
# 
# 
# kuenExtra$X.SampleID[grep("Kenya", kuenExtra$X.SampleID)]
# actualSamples[grep("Kenya", actualSamples)]
# 
# test <- c(kuenExtra2$X.SampleID, kuenExtra$X.SampleID)
# length(test)
# length(unique(test))
# length(kuenExtra$X.SampleID)
# length(kuenExtra2$X.SampleID)
# 
# 
# kuenExtra <- kuenExtra %>% filter(X.SampleID %in% actualSamples)
# Match and cross-reference location/samples so that I can add in extra metadata here


## TEMPORARILY REMOVE WALKE2015 BECAUSE NO METADATA
# toRemoveSecondaryAccess <- c("SRP062596", "SRP062876","SRP065158","SRP065432", "SRP087497", "SRP111171", "SRP115972", "SRP168079", "SRP168080", "SRP185782", "SRP268531" ,"SRP149982", "SRP168824", "SRP284492","SRP110604")

## List of final studies with title names
paperInfo <- read.delim("00_all_sample_meta/studies_to_include.txt")

## Import non-EBI
Prest <- read.delim("00_manual_metadata/PREST_meta_for_demultiplexing_FULL.txt", sep=",")
Korpita <- read.delim("00_manual_metadata/KORPITA_meta_for_demultiplexing.txt")

# Merge studies
amph_filt <- amph %>% full_join(walke_study_full) %>%
  full_join(longo_study_full) %>%
  filter(ENA.STUDY != "SRP285287") %>% # This is the ruthsatz study
  full_join(ruthsatz_study_full) #%>%
  # filter(!`ENA.STUDY` %in% toRemoveSecondaryAccess)

 
 table(amph_filt$ENA.STUDY)
allStudies <- unique(amph_filt$ENA.STUDY)
length(allStudies)

######  Go through each study and tidy up metadata######
keepCol <- c("ENASTUDY"
             , "accession"
             , "ENARUN"
             , "ENASUBMISSION"
             , "SampleID"
             , "LifeStage"
             , "SampleType"
             , "CollectionDate"
             , "Location"
             , "Bd_status"
             , "Host"
             , "Elevation"
             , "Latitude"
             , "Longitude"
             , "Captive"
             , "Pooled")
allStudies
##### SRP112703 #####
study_SRP112703 <- amph_filt %>% filter(ENA.STUDY=="SRP112703") 
study_SRP112703_filt <- study_SRP112703[, which(apply(study_SRP112703,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  mutate(Latitude = ifelse(Latitude == "not", NA, Latitude)) %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  filter(!host %in% c("Leaf", "Soil","Water","Stream.sed")) %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(one_of(keepCol)) 



##### SRP149704 #####
study_SRP149704 <- amph_filt %>% filter(ENA.STUDY=="SRP149704") 
study_SRP149704_filt <- study_SRP149704[, which(apply(study_SRP149704,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  mutate(Latitude = ifelse(Latitude == "not", NA, Latitude)) %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = NA
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(one_of(keepCol))

##### SRP075042 #####
study_SRP075042 <- amph_filt %>% filter(ENA.STUDY=="SRP075042") 
study_SRP075042_filt <- study_SRP075042[, which(apply(study_SRP075042,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(Latitude = ifelse(Latitude == "not", NA, Latitude)) %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(one_of(keepCol))

##### SRP074714 #####
study_SRP074714 <- amph_filt %>% filter(ENA.STUDY=="SRP074714") 

study_SRP074714_filt <- study_SRP074714[, which(apply(study_SRP074714,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(Elevation = ifelse(geo_loc_name == "Germany: Elm", 238, 92)) %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = Elevation
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
select(one_of(keepCol))

##### SRP078056 #####
study_SRP078056 <- amph_filt %>% filter(ENA.STUDY=="SRP078056") 

study_SRP078056_filt <- study_SRP078056[, which(apply(study_SRP078056,2,function(x) any(!is.na(x))))] %>% 
  rowwise() %>% mutate(isSkin = ifelse(length(grep("Skin", isolation_source))>0, TRUE, FALSE)) %>% ungroup() %>%
  filter(isSkin) %>%
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  mutate(Latitude = ifelse(Latitude == "not", NA, Latitude)) %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(one_of(keepCol))

##### ERP118128 #####
study_ERP118128 <- amph_filt %>% filter(ENA.STUDY=="ERP118128") 

study_ERP118128_filt <- study_ERP118128[, which(apply(study_ERP118128,2,function(x) any(!is.na(x))))] %>% 
  filter(is.na(environmental_sample), !is.na(bd_status)) %>%
  mutate(Elevation = ifelse(site == "Jones_Pond", 1219
                            , ifelse(site == "Park_Lake", 1828,
                                     ifelse(site == "Gipsy_Lake", 2243,
                                            ifelse(site == "Doney_Lake", 2170, NA ))))) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection.date
         ,Location = site
         ,Bd_status = bd_status
         ,Host = "Rana luteiventris" 
         ,Longitude =geographic.location..longitude.
          , Latitude = geographic.location..latitude.
         , Elevation = Elevation
         , Captive = FALSE
         , Pooled = FALSE
  ) %>%
  select(one_of(keepCol))

##### SRP066198 #####
study_SRP066198 <- amph_filt %>% filter(ENA.STUDY=="SRP066198") 

study_SRP066198_filt <- study_SRP066198[, which(apply(study_SRP066198,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>% 
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = host_life_stage
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = ifelse(host_disease == "Batrachochytrium_dendrobatidis","Pos","Neg")
         ,Host = host 
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Pooled=FALSE
         , Captive = TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP136053 #####
# allStudies[4]
study_SRP136053 <- amph_filt %>% filter(ENA.STUDY=="SRP136053") 

study_SRP136053_filt <- study_SRP136053[, which(apply(study_SRP136053,2,function(x) any(!is.na(x))))] %>% 
  filter(sample_type =="skin swab") %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Tadpole"
         ,SampleType = sample_type
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host 
         ,Longitude = NA
         , Latitude = NA
         , Elevation = NA
         , Pooled=FALSE
         , Captive= FALSE
  ) %>%
  select(one_of(keepCol))

##### SRP201700 #####
# NO META ON SPECIES
study_SRP201700 <- amph_filt %>% filter(ENA.STUDY=="SRP201700")

study_SRP201700_filt <- study_SRP201700[, which(apply(study_SRP201700,2,function(x) any(!is.na(x))))] %>%
  separate(alias, into=c("un1","un2","un3"), remove=FALSE, sep="-", fill="left") %>%
  # select(alias, un1, un2, un3, collection_date, isolation_source, temp) %>% View()
  filter(isolation_source=="skin", !un3 %in% c("S1","S2","S3","S4","S5")) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = "Tibetan Mountains"
         ,Bd_status = NA
         ,Host = NA
         ,Longitude = NA
         , Latitude = NA
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
    select(one_of(keepCol))

##### SRP074716 #####
# allStudies[6]
study_SRP074716 <- amph_filt %>% filter(ENA.STUDY=="SRP074716") 

study_SRP074716_filt <- study_SRP074716[, which(apply(study_SRP074716,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Larvae"
         ,SampleType = isolation_source
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled = FALSE
  ) %>%
  select(one_of(keepCol))


##### SRP130045 #####
# allStudies[7]
study_SRP130045 <- amph_filt %>% filter(ENA.STUDY=="SRP130045") 

# Lots of samples need to be removed because there was Bd exposure.
study_SRP130045_filt <- study_SRP130045[, which(apply(study_SRP130045,2,function(x) any(!is.na(x))))] %>% 
  separate(alias, into = c("sp","day"), remove=FALSE) %>%
  mutate(Captive = ifelse(day == "MB0", FALSE, TRUE)) %>%
  mutate(Bd_status = "Neg") %>% 
  filter(day %in% c("MB0","MB1")) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = paste0("UK: ",isolation_source)
         ,Bd_status = Bd_status
         ,Host = host
         ,Longitude = 0.2
         , Latitude = 52.5
         , Elevation = 250
         , Captive = Captive
         , Pooled=FALSE
  ) %>%
  select(one_of(keepCol))


##### SRP110604 #####
study_SRP110604 <- amph_filt %>% filter(ENA.STUDY=="SRP110604")

# There's an extra row here
study_SRP110604_filt <- study_SRP110604[, which(apply(study_SRP110604,2,function(x) any(!is.na(x))))] %>%
  separate(ENA.RUN, into=c("ENA.RUN", "ENA.RUN2"), sep=",") %>%
  filter(is.na(ENA.RUN2)) %>%
  mutate(Species = gsub("(_| )Cf(_| )"," ",host)) %>%
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  rowwise() %>%  mutate(year = paste0("20",substr(collection_date, 1,2)), day=substr(collection_date, 3,4), month = gsub("^.*-","",collection_date)) %>% ungroup() %>%
  unite(year, month, day, col=collection_date, sep = "-") %>% 
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = isolation_source
         ,Bd_status = NA
         ,Host = "Rhinella marina"
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))

##### SRP133360 #####
### NOTE: MIXED INFECTION STATUS here too
study_SRP133360 <- amph_filt %>% filter(ENA.STUDY=="SRP133360") 

study_SRP133360_filt <- study_SRP133360[, which(apply(study_SRP133360,2,function(x) any(!is.na(x))))] %>% 
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Metamorph"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = "Acris blanchardi"
         ,Longitude = NA
         , Latitude = NA
         , Elevation = NA
         , Captive = TRUE
         , Pooled = TRUE
  ) %>%
  select(all_of(keepCol))

#### SRP097778 #####
# BOTH LARVAL AND ADULTS INCLUDED
study_SRP097778 <- amph_filt %>% filter(ENA.STUDY=="SRP097778") 

study_SRP097778_filt <- study_SRP097778[, which(apply(study_SRP097778,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = NA
          ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(all_of(keepCol))

##### SRP253593 REMOVED #####
## ARE THESE HISEQ? Yes. REMOVE
# study_SRP253593 <- amph_filt %>% filter(ENA.STUDY=="SRP253593") 
# 
# study_SRP253593_filt <- study_SRP253593[, which(apply(study_SRP253593,2,function(x) any(!is.na(x))))] %>% 
#   filter(env_broad_scale=="Skin" | env_broad_scale=="Stomach content") %>%
#   mutate(ENASTUDY = ENA.STUDY
#          ,ENARUN=ENA.RUN
#          , ENASUBMISSION=ENA.SUBMISSION
#          , SampleID = alias
#          ,LifeStage = "Adult"
#          ,SampleType = ifelse(env_broad_scale=="Skin", "skin swab", "Stomach content")
#          ,CollectionDate = collection_date
#          ,Location = geo_loc_name
#          ,Bd_status = NA
#          ,Host = host
#          ,Longitude = NA
#          , Latitude = NA
#          , Elevation = NA
#          , Captive = NA
#          , Pooled = FALSE
#   ) %>%
#   select(all_of(keepCol))


##### SRP151078 #####
# allStudies[14]
study_SRP151078 <- amph_filt %>% filter(ENA.STUDY=="SRP151078") 

elevation_manual <- data.frame(geo_loc_name = c("Germany: Solling"
                                            , "Germany: Eifel, Zweifallshammer"
                                            , "Germany: Eifel,Fleissbach"
                                            , "Germany: Eifel, Kallerbach"
                                            , "Germany: Kottenforst"
                                            , "Germany: Harz"
                                            , "Germany: Eifel, Rosbach"
                                            , "Belgium: Ughent Lab"
                                            ,"Germany: Eifel, Weisse Wehe"
                                            ,"Germany: Eifel, Sauerbach"
                                            ,"Germany: Eifel, Haftenbach"
                                            ,"Germany: Eifel, Solchbachtal"
                                            ,"Germany: Eifel,ER"
                                            ,"Germany: Rhineland-Palatinate, 6104")
                               , Elevation = c(360, 240, 260, 330, 150, 310, 275, 357, 330, 410, 400, 330, 270, 320))

study_SRP151078_filt <- study_SRP151078[, which(apply(study_SRP151078,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  left_join(elevation_manual) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = "Neg"
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = Elevation
         , Captive = FALSE
         , Pooled=FALSE
  ) %>%
  select(all_of(keepCol))


##### SRP075048 #####
# allStudies[15]
study_SRP075048 <- amph_filt %>% filter(ENA.STUDY=="SRP075048") 

study_SRP075048_filt <- study_SRP075048[, which(apply(study_SRP075048,2,function(x) any(!is.na(x))))] %>% 
  mutate(lat_lon = ifelse(lat_lon=="Not applicable",NA, lat_lon)) %>%
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = isolation_source
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = "Neg"
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = TRUE
         , Pooled = FALSE
  ) %>%
  select(all_of(keepCol))

##### SRP065432 #####
## POOLED ###
study_SRP065432 <- amph_filt %>% filter(ENA.STUDY=="SRP065432") 

study_SRP065432_filt <- study_SRP065432[, which(apply(study_SRP065432,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         , ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         , LifeStage = NA
         ,SampleType = "skin swab"
         ,CollectionDate = NA
         ,Location = NA
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive = FALSE
         , Pooled=TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP245448 #####
# allStudies[16]
study_SRP245448 <- amph_filt %>% filter(ENA.STUDY=="SRP245448") 

study_SRP245448_filt <- study_SRP245448[, which(apply(study_SRP245448,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  filter(is.na(isolation_source)) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = isolation_source
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Pooled=FALSE
         , Captive= NA
  ) %>%
  select(all_of(keepCol))

##### SRP185782 #####
# FORWARD REVERSE; Can't demultiplex
study_SRP185782 <- amph_filt %>% filter(ENA.STUDY=="SRP185782") 

study_SRP185782_filt <- study_SRP185782[, which(apply(study_SRP185782,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = NA
         ,SampleType = NA
         ,CollectionDate = NA
         ,Location = NA
         ,Bd_status = NA
         ,Host = NA
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Pooled=TRUE
         , Captive= NA
  ) %>%
  select(all_of(keepCol))


##### SRP242301 #####
### WORKING HERE 
# SAMPLES WERE POOLED; NEED TO SEPARATE OUT
study_SRP242301 <- amph_filt %>% filter(ENA.STUDY=="SRP242301")

allAccTemp <- study_SRP242301 %>% pull(ENA.RUN) %>% unique()
allRunAcc <- study_SRP242301 %>% select(accession,EXTERNAL_ID, ENA.RUN) %>% 
  rename(ENA_grouped = ENA.RUN) %>%
  distinct()
for ( a in allAccTemp) {
  # a=allAccTemp[1]
  subGroup <- unlist(lapply(a, FUN=function(x) strsplit(x, split=",")))
  
  listSubGroup <- c()
  for ( b in subGroup) {
    if (length(grep("-", b))>0) {
      ends <- unlist(strsplit(b, split="-"))
      start_unlisted <- unlist(strsplit(ends[1],""))
      end_unlisted <- unlist(strsplit(ends[2],""))
      digitsVary <- which(start_unlisted != end_unlisted)
      constantString <- substr(ends[1], 1,digitsVary[1]-1)
      endings <- seq(as.numeric(substr(ends[1], digitsVary[1], nchar(ends[1]))),as.numeric(substr(ends[2], digitsVary[1], nchar(ends[2]))))
      listSubGroup <- c(listSubGroup, paste0(constantString, endings))
      
    } else {
      listSubGroup <- c(listSubGroup, b)
    }
  }
  allRunAcc <- full_join(allRunAcc, data.frame(ENA_grouped = a, ENA.RUN = listSubGroup))
}

study_SRP242301_filt <- study_SRP242301[, which(apply(study_SRP242301,2,function(x) any(!is.na(x))))] %>%
  rename(ENA_grouped = ENA.RUN) %>%
  full_join(allRunAcc) %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  filter(host!="water") %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION = ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA 
         , Captive = FALSE
         , Pooled = FALSE
  ) %>%
  select(all_of(keepCol))


##### SRP149982 #####
study_SRP149982 <- amph_filt %>% filter(ENA.STUDY=="SRP149982")

study_SRP149982_filt <- study_SRP149982[, which(apply(study_SRP149982,2,function(x) any(!is.na(x))))] %>%
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  filter(!is.na(host)) %>%
  mutate(altitude = as.numeric(gsub(" m","",altitude))) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         ,LifeStage = host_life_stage
         ,SampleType = ifelse(env_material == "skin epidermis", "skin swab", env_material)
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = altitude
         , Captive = FALSE
         , Pooled = FALSE
  ) %>%
select(keepCol)

##### SRP104695 #####
# allStudies[19]
study_SRP104695 <- amph_filt %>% filter(ENA.STUDY=="SRP104695") 

study_SRP104695_filt <- study_SRP104695[, which(apply(study_SRP104695,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(host = gsub("_"," ", host))%>%
  filter(!is.na(host)) %>%
  rowwise() %>% mutate(ComboDate=paste( collection_date, Date, sep="-")) %>% ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = X.SampleID
         ,LifeStage = ifelse(Ageclass == "Ad", "Adult", ifelse(Ageclass=="juv_sub","Subadult",Ageclass))
         ,SampleType = "skin swab"
         ,CollectionDate = ComboDate
         ,Location = paste0(geo_loc_name, Population)
         ,Bd_status = ifelse(Bd==1, "Pos", ifelse(Bd==0,"Neg", NA))
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled = TRUE
  ) %>%
select(one_of(keepCol))


##### SRP087497 #####
study_SRP087497 <- amph_filt %>% filter(ENA.STUDY=="SRP087497") 

study_SRP087497_filt <- study_SRP087497[, which(apply(study_SRP087497,2,function(x) any(!is.na(x))))] %>% 
 filter(!is.na(host)) %>%
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         , LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = host
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled = TRUE
  ) %>%
select(one_of(keepCol))


##### SRP065158 #####
# allStudies[20]
study_SRP065158 <- amph_filt %>% filter(ENA.STUDY=="SRP065158") 
study_SRP065158_filt <- study_SRP065158[, which(apply(study_SRP065158,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         , LifeStage = NA
         ,SampleType = NA
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         ,Host = NA
         ,Longitude = Longitude
         , Latitude = Latitude
         , Elevation = NA
         , Captive = FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))



##### SRP062395 #####
study_SRP062395 <- amph_filt %>% filter(ENA.STUDY=="SRP062395") 

study_SRP062395_filt <- study_SRP062395[, which(apply(study_SRP062395,2,function(x) any(!is.na(x))))] %>% 
  filter(!is.na(host)) %>%
  mutate(host_scientific_name = ifelse(host=="Bullfrog", "Lithobates actesbeianus"
                                       , ifelse(host=="Newt","Notophthalmus viridescens"
                                                , ifelse(host=="Peeper", "Pseudacris crucifer"
                                                         , ifelse(host=="American Toad", "Anaxyrus americanus", NA))))) %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = alias
         , LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_timestamp
         ,Location = paste0(country, ":", site)
         ,Bd_status = NA
         ,Host = host_scientific_name
         ,Longitude = longitude
         , Latitude = latitude
         , Elevation = elevation
         , Captive = FALSE
         , Pooled = TRUE
  ) %>%
select(all_of(keepCol))



##### SRP111171 #####
study_SRP111171 <- amph_filt %>% filter(ENA.STUDY=="SRP111171") 

study_SRP111171_filt <- study_SRP111171[, which(apply(study_SRP111171,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))




##### SRP062876 #####
study_SRP062876 <- amph_filt %>% filter(ENA.STUDY=="SRP062876") 
study_SRP062876_filt <- study_SRP062876[, which(apply(study_SRP062876,2,function(x) any(!is.na(x))))] %>% 
    separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = "Adult"
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP062596 #####
study_SRP062596 <- amph_filt %>% filter(ENA.STUDY=="SRP062596") 
study_SRP062596_filt <- study_SRP062596[, which(apply(study_SRP062596,2,function(x) any(!is.na(x))))] %>% 
    mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = NA
         ,SampleType =NA
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = NA
         , Latitude = NA
         ,Host = NA 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))



##### SRP168080 #####
study_SRP168080 <- amph_filt %>% filter(ENA.STUDY=="SRP168080") 
study_SRP168080_filt <- study_SRP168080[, which(apply(study_SRP168080,2,function(x) any(!is.na(x))))] %>% 
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = NA
         ,SampleType =NA
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = NA
         , Latitude = NA
         ,Host = NA 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))

##### SRP168079 #####
study_SRP168079 <- amph_filt %>% filter(ENA.STUDY=="SRP168079") 
study_SRP168079_filt <- study_SRP168079[, which(apply(study_SRP168079,2,function(x) any(!is.na(x))))] %>% 
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = NA
         ,SampleType =NA
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = NA
         , Latitude = NA
         ,Host = NA 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP268531 #####
study_SRP268531 <- amph_filt %>% filter(ENA.STUDY=="SRP268531") 
study_SRP268531_filt <- study_SRP268531[, which(apply(study_SRP268531,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = NA
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP115972 #####
study_SRP115972 <- amph_filt %>% filter(ENA.STUDY=="SRP115972") 
study_SRP115972_filt <- study_SRP115972[, which(apply(study_SRP115972,2,function(x) any(!is.na(x))))] %>% 
  separate(lat_lon, into = c("Latitude","d1","Longitude","d2"), sep = " ") %>%
  rowwise() %>%
  mutate(Latitude=ifelse(d1=="S", -1*as.numeric(Latitude), as.numeric(Latitude)), Longitude=ifelse(d2=="W", -1*as.numeric(Longitude), as.numeric(Longitude))) %>%
  ungroup() %>%
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = NA
         ,SampleType = "skin swab"
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = NA
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = TRUE
  ) %>%
  select(one_of(keepCol))


##### SRP285287 #####
study_SRP285287 <- amph_filt %>% filter(ENA.STUDY=="SRP285287") 
long_manual <-data.frame(geo_loc_name = c("Brazil: Rio de Janeiro, Teresopolis"
                                          ,"Brazil: Parana, Morretes"
                                          ,"Brazil: Santa Catarina, Sao Miguel d Oeste"
                                          ,"Brazil: Santa Catarina, Pomerode"
                                          ,"Brazil: Espirito Santo, Linhares"
                                          ,"Brazil: Sao Paulo, Itariri"
                                          ,"Brazil: Sao Paulo, Iporanga"
                                          ,"Brazil: Bahia, Mata de Sao Joao"
                                          ,"Brazil: Espirito Santo, Vargem Alta"
                                          
                                          , "Brazil: Sao Paulo, Ubatuba"
                                          ,"Brazil: Alagoas, Murici"
                                          , "Brazil: Espirito Santo, Santa Teresa"
                                          , "Brazil: Sao Paulo, Jundiai"
                                          ,"Brazil: Bahia, Camacan"
                                          ,"Brazil: Rio Grande do Norte, Baia Formosa"
                                          ,"Brazil: Sao Paulo, Iguape"
                                          ,"Brazil: Santa Catarina, Rancho Queimado"
                                          ,"Brazil: Paraiba, Joao Pessoa"
                                          , "Brazil: Sao Paulo, Pedro de Toledo"
                                          ,"Brazil: Rio Grande do Sul, Torres"
                                          , "Brazil: Sergipe, Areia Branca"
                                          ,"Brazil: Pernambuco, Paulista"
                                          ,"Brazil: Rio de Janeiro, Macae, Barra do Sana")
                         , Longitude = c(-42.9, -48.8, -53.5, -49.1, -40, -47.1, -48.5, -38.3, -41
                                         , -45, -35.9, -40.6, -46, -39.4, -35, -47.5, -49, -34.8, -47.2, -49.8, -37.3, -34.9, -42.2))

study_SRP285287_filt <- study_SRP285287[, which(apply(study_SRP285287,2,function(x) any(!is.na(x))))] %>% 
  left_join(long_manual) %>% 
  # filter(host=="Dendropsophus elegans") %>% select(geo_loc_name, State, Municipality, Group) %>% distinct() %>%View()
  mutate(ENASTUDY = ENA.STUDY
         ,ENARUN=ENA.RUN
         , ENASUBMISSION=ENA.SUBMISSION
         , SampleID = SUBMITTER_ID.label.Sample.name
         ,LifeStage = dev_stage
         ,SampleType = tissue
         ,CollectionDate = collection_date
         ,Location = geo_loc_name
         ,Bd_status = ifelse(Bd.status=="+","Pos",ifelse(Bd.status=="-","Neg",NA))
         , Longitude = Longitude
         , Latitude = Latitude
         ,Host = host 
         , Elevation = Elevation
         , Captive=FALSE
         , Pooled = FALSE
  ) %>%
  select(one_of(keepCol))

#### Prest ####
study_Prest_filt <- Prest%>% 
  filter(Species !="Environment", Sample_Type !="Skin side #2", Other_Notes !="Dead") %>%
  rowwise() %>% mutate(LifeStage  = ifelse(length(grep("Adult", Lifestage))>0, "Adult", 
                                          ifelse(length(grep("(m|M)eta", Lifestage))>0, "Metamorph",
                                                 ifelse(length(grep("(t|T)adpole", Lifestage))>0, "Tadpole", 
                                                        ifelse(length(grep("(j|J|Je)uv", Lifestage))>0, "Juvenile", Lifestage))))) %>% ungroup() %>% 
  mutate(Species = gsub("_", " ", Species)) %>%
  rowwise() %>% mutate(SampleType = ifelse(length(grep("Body", Sample_Type))>0, "skin swab",
                                           ifelse(length(grep("(w|W)ater", Sample_Type))>0, "Environ",
                                                  ifelse(length(grep("(s|S)kin", Sample_Type))>0, "skin swab",
                                                         ifelse(Sample_Type=="", "skin swab",Sample_Type))))) %>% ungroup() %>%
  mutate(Location = paste0("USA:",ifelse(Location=="4 mile","4 Mile", Location))) %>%
  mutate(ENASTUDY = "Prest"
         ,ENARUN=Sample.id
         , ENASUBMISSION=NA
         , accession = Sample.id
         , SampleID = Sample.id
         ,LifeStage = LifeStage 
         ,SampleType = SampleType
         ,CollectionDate = Date
         ,Location = Location
         ,Bd_status = NA
         , Longitude = NA
         , Latitude = NA
         ,Host = Species 
         , Elevation = NA
         , Captive=FALSE
         , Pooled = FALSE
  ) %>%
select(one_of(keepCol))

#### Korpita #####
study_Korpita_filt <- Korpita%>% 
  filter(SampleType=="Swab", Species == "BToad") %>%
  mutate(Host = ifelse(Species=="BToad", "Anaxyrus boreas", Species) ) %>%
  mutate(SampleType = ifelse(BodySite %in% c("Skin","Swab"), "skin swab", BodySite)) %>%
  rowwise() %>%
  mutate(nonTreat = ifelse(length(grep("Pre",Treatment))>0, "keep", 
                             ifelse(is.na(Treatment), "keep","not"))) %>%
  ungroup() %>% 
  filter(nonTreat == "keep") %>%
  mutate(LocationCollected = paste0("USA:",LocationCollected)) %>%
  mutate(ENASTUDY = "Korpita"
         ,ENARUN=X.SampleID
         , ENASUBMISSION=NA
         , accession = X.SampleID
         , SampleID = X.SampleID
         ,LifeStage = LifeStage 
         ,SampleType = SampleType
         ,CollectionDate = DateCollected
         ,Location = LocationCollected
         ,Bd_status = NA
         , Longitude = NA
         , Latitude = NA
         , Host = Host 
         , Elevation = NA
         , Captive=ifelse(CaptiveWild=="Captive", TRUE, FALSE)
         , Pooled = FALSE
  ) %>%
  select(one_of(keepCol))


#######################################

#### Combine all individual datasets #####
all_meta <- as.data.frame(matrix(ncol=length(keepCol), nrow=0, dimnames = list(NULL, keepCol)))
for ( s in c(allStudies, "Prest", "Korpita")) {
  all_meta <- rbind(all_meta, get(paste0("study_",s,"_filt")))
}

unique(all_meta$ENASTUDY)

##### Cross-reference with Jordan's samples ##### 
#FUll sample
allKuen <- read.delim("00_nonEBI_data/Kuen/summary_og_table.txt", sep=":")
allKuenSamp <- allKuen$Num.samples[13:nrow(allKuen)]
kuenFilt1 <- kuenExtra %>% filter(X.SampleID %in% allKuenSamp)
kuenFilt2 <- kuenExtra %>% filter(!X.SampleID %in% allKuenSamp) %>%
  rowwise() %>% filter(length(grep("MVB", X.SampleID))>0) %>% 
  mutate(X.SampleID = gsub("MVB-", "MVB", X.SampleID)) %>%
  filter(X.SampleID %in% allKuenSamp)
kuen_filtered_meta <- rbind(kuenFilt1,kuenFilt2)
  
kuen_compare <- kuen_filtered_meta %>% select(-FrogSpecies.1, -LifeStage.1) %>%
  mutate(Captive.1 = ifelse(Captive_Wild=="Captive", TRUE, FALSE)) %>%
  mutate(FrogSpecies = ifelse(FrogSpecies == "Boophis laurenti ","Boophis laurenti", FrogSpecies)) %>%
  mutate(FrogSpecies = ifelse(FrogSpecies == "Guibemantis timidus ", "Guibemantis timidus", FrogSpecies)) %>%
  mutate(FrogSpecies = ifelse(FrogSpecies == "mantella laevigata", "Mantella laevigata", FrogSpecies)) %>%
  mutate(SampleType.1 = ifelse(SampleType%in% c("FullBody","Dorsal","Ventral"), "skin swab", SampleType)) %>%
  filter(SampleType.1 =="skin swab") %>%
  unite(SequenceCenter,SampleContact, col="Study") %>%
  rowwise() %>%
  mutate(Location.1=paste0(Country, ":", Locality)) %>% ungroup() %>%
  mutate(Location.1 = ifelse(Location.1=="USA:Sawhill ponds","USA:Sawhill",Location.1)) %>%
  mutate(KUEN=TRUE) %>%
  rename(LifeStage.1 = LifeStage, SampleID = X.SampleID, Latitude.1=Latitude, Longitude.1=Longitude, Host.1 = FrogSpecies, CollectionDate.1 = CollectionDate, Elevation.1=Revised_Altitude) %>%
  select(KUEN, SampleID, Host.1, Latitude.1, Longitude.1, LifeStage.1, Location.1, Captive.1, CollectionDate.1, Elevation.1, SampleType.1, Study) %>%
  full_join(all_meta %>% mutate(Study=ENASTUDY)) %>% mutate(KUEN=ifelse(is.na(KUEN), FALSE, KUEN), Pooled=ifelse(is.na(Pooled), FALSE, Pooled)) %>%  
  mutate(LifeStage = ifelse(KUEN & (is.na(LifeStage) | LifeStage!=LifeStage.1) & !is.na(LifeStage.1), LifeStage.1, LifeStage)) %>%
  mutate(Latitude = ifelse(KUEN & (is.na(Latitude)) & !is.na(Latitude.1), Latitude.1, Latitude)) %>%
  mutate(Longitude = ifelse(KUEN & (is.na(Longitude)) & !is.na(Longitude.1), Longitude.1, Longitude))  %>%
  mutate(Elevation = ifelse(KUEN & (is.na(Elevation)) & !is.na(Elevation.1), Elevation.1, Elevation))  %>%
  mutate(Host = ifelse(KUEN & (is.na(Host)) & !is.na(Host.1), Host.1, Host))  %>%
  mutate(Location = ifelse(KUEN & (is.na(Location)) & !is.na(Location.1), Location.1, Location))  %>%
  # mutate(Location = ifelse(Location=="USA:Sawhillponds", "USA:Sawhill", Location)) %>%
  mutate(CollectionDate = ifelse(KUEN & (is.na(CollectionDate)) & !is.na(CollectionDate.1), CollectionDate.1, CollectionDate)) %>%
  mutate(SampleType = ifelse(KUEN & (is.na(SampleType)) & !is.na(SampleType.1), SampleType.1, SampleType))

new_all_meta <- kuen_compare %>%
  select(one_of(keepCol), KUEN, Study) %>%
   mutate(ENASTUDY = ifelse(KUEN & is.na(ENASTUDY), "Kuen", ENASTUDY))
unique(new_all_meta$Study)
# new_all_meta %>% filter(ENASTUDY=="Kuen")
# 
# askJordan <- kuen_compare[,order(colnames(kuen_compare))] %>% filter(KUEN) %>%
#   mutate(sameHost = Host==Host.1) %>%
#   select(CollectionDate, Location, sameHost, Host, Host.1) %>%
#   mutate(sameHost=ifelse(is.na(sameHost), FALSE, sameHost)) %>%
#   filter(Host%in% c("Boophis quasiboehmei", "Mantella baroni", "Mantidactylus curtus", "Mantidactylus femoralis", "Mantidactylus grandidieri", "Salamandra salamandra", "Gephyromantis luteus")) %>%
#   filter(!sameHost) %>% 
#   distinct %>% arrange(sameHost, Host) %>% select(CollectionDate, Location, Host, Host.1)
# 
# write.table(askJordan, file="01_parsed_meta/askJordan.txt", quote=FALSE, row.names = FALSE, sep="\t")

## Filter out pooled samples 
new_all_meta %>% filter(!Pooled) %>% pull(Study) %>% unique()

# Add in latlong of shared sites between mckenzie lab studies-- some manual spots too
# manual spots 
new_all_meta %>%
  mutate(Location=gsub(" ","",Location)) %>%
  filter(Study %in% c("Colorado_Jordan.Kueneman","Prest","Korpita"))%>%
  # filter(length(grep("wood", Location))>0) %>%
  select(Location, Latitude, Longitude) %>% distinct() %>% filter(is.na(Latitude)) %>%
  distinct()
newDat <- rbind(c("USA:BoulderCreek", 40.31923835373663, -105.60685923580841)
      , c("USA:NASRF", 37.47736227284562, -105.94477652064859)
        , c("USA:LostLake", 40.50830310618662, -105.60012066546804)
      , c("USA:Spruce", 40.341899091434826, -105.68679131148679)
      , c("USA:Zimmerman", 40.53872994427544, -105.88322456395508)
      , c("USA:4Mile", 38.974698135367376, -106.14546846521489)) %>% as.data.frame() %>%
  rename(Location=V1, Latitude=V2, Longitude=V3) %>%
  mutate(Latitude=as.numeric(Latitude), Longitude=as.numeric(Longitude))


mckenzie_coords <- new_all_meta %>%
  mutate(Location=gsub(" ","",Location)) %>%
  filter(Study %in% c("Colorado_Jordan.Kueneman","Prest","Korpita"))%>%
  # filter(length(grep("wood", Location))>0) %>%
  select(Location, Latitude, Longitude) %>% distinct() %>%
  drop_na() %>% mutate(Latitude=as.numeric(Latitude), Longitude=as.numeric(Longitude)) %>%
  full_join(newDat)
  
# all_meta %>% filter(Host=="Anaxyrus boreas") %>% View()
all_meta_filt <- new_all_meta %>%
  mutate(Location=gsub(" ","",Location), Longitude=str_trim(Longitude), Latitude=str_trim(Latitude)) %>%
  mutate(Longitude=as.numeric(Longitude), Latitude=as.numeric(Latitude)) %>%
  rowwise() %>%
  mutate(Latitude=ifelse(Location %in% mckenzie_coords$Location, pull(mckenzie_coords[match(Location, mckenzie_coords$Location),"Latitude"]), Latitude)
         ,Longitude=ifelse(Location %in% mckenzie_coords$Location, pull(mckenzie_coords[match(Location, mckenzie_coords$Location),"Longitude"]), Longitude)) %>% 
  ungroup() %>%
  filter(!is.na(ENASTUDY), !Pooled, LifeStage == "Adult", SampleType == "skin swab", !is.na(Host)) %>%
  rowwise() %>% mutate(unknownsp = length(grep("(sp|spp)[.]$",Host))>0) %>% ungroup() %>%
  rowwise() %>% mutate(cfsp = length(grep("cf.",Host))>0) %>% ungroup() %>%
  filter(!unknownsp, !cfsp , !Host %in% c("Pseudoeurycea sp. n. Mozotal 2"
                                  ,"Bolitoglossa franklini lincolni"
                                  ,"Craugastor"
                                  ,"Pseudoeurycea sp. n. Mozotal 1"
                                  ,"Scinax aromothylla"
                                  ,"Scinax gr. ruber")) %>%
  mutate(Host = gsub(" x ", " ", Host)) %>%
  mutate(Host = gsub("Agalychnis granulosa", "Hylomantis granulosa", Host)) %>%
  mutate(Host = gsub("Aplastodiscus eugenoi", "Aplastodiscus eugenioi", Host)) %>%
  mutate(Host = gsub("Aplastodiscus pervidis", "Aplastodiscus perviridis", Host)) %>%
  mutate(Host = gsub("Bolitoglossa franklini nigroflavescens", "Bolitoglossa franklini", Host)) %>%
  mutate(Host = gsub("Lithobates actesbeianus", "Rana catesbeiana", Host)) %>%
  mutate(Host = gsub("Lithobates catesbeianus", "Rana catesbeiana", Host)) %>%
  mutate(Host = gsub("Lithobates catesbeiana", "Rana catesbeiana", Host)) %>%
  mutate(Host = gsub("Lithobates clamitans", "Rana clamitans", Host)) %>%
  mutate(Host = gsub("Lithobates sylvaticus", "Rana sylvatica", Host)) %>%
  mutate(Host = gsub("Pithecopus rhodei", "Phyllomedusa rohdei", Host)) %>%
  mutate(Host = gsub("Scinax flavogutattus", "Scinax flavoguttatus", Host)) %>%
  mutate(Host = gsub("Scinax flavogutattus", "Scinax flavoguttatus", Host)) %>%
  mutate(Host = gsub("Boophis goudoti", "Boophis goudotii", Host)) %>%
  mutate(Host = gsub("Anodonthyla boulengeri", "Anodonthyla boulengerii", Host)) %>%
  mutate(Host = gsub("Cophyla grandis", "Platypelis grandis", Host)) %>%
  mutate(Host = gsub("Boophis boehemi", "Boophis boehmei", Host)) %>%
  mutate(Host = gsub("R. cat", "Rana catesbeiana", Host)) %>%
  mutate(Host = gsub("Mantidactylus cowanii small", "Mantidactylus cowanii", Host)) %>%
  mutate(Host = gsub("Hylomantis lemur", "Agalychnis lemur", Host)) %>%
  mutate(Host = gsub("Eleutherodactylus caspari", "Eleutherodactylus casparii", Host)) %>%
  mutate(Host = gsub("iii", "ii", Host)) # IN case there was added extra
  

rough_table <- all_meta_filt %>% 
  group_by(Host) %>% summarise(nSamples=n(), nStudies = length(unique(ENASTUDY)), nLocations = length(unique(Location)), nTimepoints = length(unique(CollectionDate))) %>% 
  arrange(-nSamples, -nStudies, -nLocations)
nrow(rough_table)
nrow(all_meta_filt)
all_meta_filt %>% filter(Study %in% c("Kuen","Prest","Korpita")) %>% select(Location, Latitude, Longitude) %>% distinct()

write.table(rough_table, quote=FALSE, sep="\t", row.names = FALSE, file = "01_parsed_meta/rough_table_of_samples.txt")

#### Get list of samples to download off of EBI/ENA
# ERR_toDownload <- all_meta_filt %>%
#   select(ENARUN, ENASUBMISSION) %>% distinct() %>%
#   rowwise() %>%
#   filter(length(grep("ERR", ENARUN))>0) %>% ungroup() %>%
#   select(ENARUN) %>% pull(ENARUN)
# 
# # ERR_folder <- unique(sapply(ERR_toDownload, function(x) substr(x,1,6)))
# 
# SRR_toDownload <- all_meta_filt %>%
#   select(ENARUN, ENASUBMISSION) %>% distinct() %>% 
#   rowwise() %>%
#   filter(length(grep("SRR", ENARUN))>0) %>% ungroup() %>%
#   select(ENARUN) %>% pull(ENARUN)
# 
# write.table(ERR_toDownload, file = "01_parsed_meta/ERR_toDownload_sample_accessions.txt", row.names = FALSE, quote = FALSE, sep="\t")
# write.table(SRR_toDownload, file = "01_parsed_meta/SRR_toDownload_sample_accessions.txt", row.names = FALSE, quote = FALSE, sep="\t")

### Filtered version
multipleStudyHosts <- all_meta_filt %>% 
  group_by(Host) %>% summarise(nSamples=n(), nStudies = length(unique(ENASTUDY)), nLocations = length(unique(Location)), nTimepoints = length(unique(CollectionDate))) %>% 
  arrange(-nSamples, -nStudies, -nLocations) %>% 
  filter(nLocations>2 | nStudies>1) %>% pull(Host)

# all_meta_filt %>% filter(Host == "Anaxyrus boreas") %>% View()
multipleStudyHosts_COMPLETE <- all_meta_filt %>% 
  unite(Host, Location, CollectionDate, col="UniqueHostLocTime", remove=FALSE) %>%
  group_by(UniqueHostLocTime) %>% mutate(nUnique=n()) %>% ungroup() %>%
  group_by(Host) %>% summarise(nSamples=n(), nStudies = length(unique(ENASTUDY)), nLocations = length(unique(Location)), nTimepoints = length(unique(CollectionDate)), nUnique=max(nUnique)) %>% ungroup() %>%
  arrange(-nSamples, -nStudies, -nLocations) %>%  
  filter(nLocations>2 | nStudies>2 | nUnique>10 ) %>% pull(Host)

sort(multipleStudyHosts_COMPLETE)
# 
# ERR_toDownload <- all_meta_filt %>% 
#   filter(Host %in% c(multipleStudyHosts), LifeStage=="Adult", SampleType=="skin swab") %>%
#   select(ENARUN, ENASUBMISSION) %>% distinct() %>%
#   rowwise() %>%
#   filter(length(grep("ERR", ENARUN))>0) %>% ungroup() %>%
#   select(ENARUN) %>% pull(ENARUN)

# ERR_folder <- unique(sapply(ERR_toDownload, function(x) substr(x,1,6)))

SRR_toDownload <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts), LifeStage=="Adult", SampleType=="skin swab") %>%
  select(ENARUN, ENASUBMISSION) %>% distinct() %>% 
  filter(ENARUN !="SRR5643035") %>%
  rowwise() %>%
  filter(length(grep("RR", ENARUN))>0) %>% ungroup() %>%
  select(ENARUN) %>% pull(ENARUN)
# View(data.frame(SRR_toDownload))
# write.table(ERR_toDownload, file = "01_parsed_meta/ERR_toDownload_sample_accessions.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(SRR_toDownload, file = "01_parsed_meta/SRR_toDownload_sample_accessions.txt", row.names = FALSE, quote = FALSE, sep="\t", col.names = FALSE)


SRR_toDownload_COMPLETE <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts_COMPLETE), LifeStage=="Adult", SampleType=="skin swab") %>%
  select(ENARUN, ENASUBMISSION) %>% distinct() %>% 
  filter(ENARUN !="SRR5643035") %>%
  rowwise() %>%
  filter(length(grep("RR", ENARUN))>0) %>% ungroup() %>%
  select(ENARUN) %>% pull(ENARUN)
nonEBI_toDownload_COMPLETE <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts_COMPLETE), LifeStage=="Adult", SampleType=="skin swab") %>%
  select(ENASTUDY, ENARUN) %>% distinct() %>% 
  filter(ENARUN !="SRR5643035") %>%
  rowwise() %>%
  filter(length(grep("RR", ENARUN))==0) %>% ungroup() %>%
  select(ENASTUDY, ENARUN)
Kuen_toAdd_COMPLETE <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts_COMPLETE), LifeStage=="Adult", SampleType=="skin swab") %>%
  filter(ENASTUDY=="Kuen") %>%
  select(SampleID)
# View(data.frame(nonEBI_toDownload_COMPLETE))
# write.table(ERR_toDownload, file = "01_parsed_meta/ERR_toDownload_sample_accessions.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(SRR_toDownload_COMPLETE, file = "01_parsed_meta/SRR_toDownload_sample_accessions_COMPLETE.txt", row.names = FALSE, quote = FALSE, sep="\t", col.names = FALSE)
write.table(nonEBI_toDownload_COMPLETE, file = "01_parsed_meta/nonEBI_toDownload_sample_accessions_COMPLETE.txt", row.names = FALSE, quote = FALSE, sep="\t", col.names = FALSE)
write.table(Kuen_toAdd_COMPLETE, file = "01_parsed_meta/Kuen_toAdd_COMPLETE.txt", row.names = FALSE, quote = FALSE, sep="\t", col.names = TRUE)

########## Save metadata
test_meta <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts), LifeStage=="Adult", SampleType=="skin swab") 

write.table(test_meta, file="01_parsed_meta/test_meta.txt", quote=FALSE, row.names = FALSE, sep="\t")

## Get all sample names
# allFileNames <- gsub("-table.qza","",list.files("02_download_ENAEBI_data_amphibian/qiime_output/tables/"))

full_meta <- all_meta_filt %>%
  filter(Host %in% c(multipleStudyHosts_COMPLETE), LifeStage=="Adult", SampleType=="skin swab") %>%
  rename(sampleid_descrip=SampleID) %>%
  mutate(sampleid=ifelse(is.na(ENARUN), sampleid_descrip, ENARUN)) %>%
  # filter(sampleid %in% allFileNames) %>%
  select(sampleid, everything())
write.table(full_meta, file="01_parsed_meta/full_meta.txt", quote=FALSE, row.names = FALSE, sep="\t")

sink("01_parsed_meta/LOG.txt")
print("Number of studies")
length(unique(full_meta$Study))
print("Number of amphibians")
length(unique(full_meta$Host))
sink()
## Write meta with proper names?
allSamples <- full_meta$sampleid
allSamples_nounscore <- gsub("_","-", allSamples)
write.table(allSamples_nounscore, file="01_parsed_meta/list_samples_to_keep_COMPLETE.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# ######## Edit metadata to add proper missing info
# unique(test_meta$ENASTUDY) %in% paperInfo$SecondaryAccession
# View(test_meta)
# test_meta %>% rename(SecondaryAccession = ENASTUDY) %>%
#   left_join(paperInfo)%>%
#   select(PrimaryAccession, SecondaryAccession, Title, FirstAuthor, Year, Journal) %>%
#   distinct() %>% View( title = "StudySummary")
# ### ADD MANUAL INFO
# test_meta %>% select(ENASTUDY, Location, Latitude, Longitude, Captive, Bd_status, Elevation ) %>%
#   distinct() %>% View(title = "KnownSoFar")


### SANITY TEST ###
# TOTAL samples
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()

# ALL Jordan samples
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()
# ALL non-Jordan samples
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(!KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()

# Samples that already exist in dataset
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(!is.na(ENARUN), KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()
# kuen_compare %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(!is.na(ENARUN), KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(ENASTUDY !="Kuen", KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()

# Samples added
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(ENASTUDY =="Kuen", LifeStage=="Adult", SampleType=="skin swab") %>% nrow()
all_meta_filt %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(ENASTUDY =="Kuen",KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()

# Number of Kuen samples added
nrow(Kuen_toAdd_COMPLETE)
# Show all samples are within dataset
any(!Kuen_toAdd_COMPLETE$SampleID %in% new_all_meta$SampleID)
# Identify samples that are different between Kuen and new all meta

new_all_meta %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(is.na(ENARUN), KUEN, LifeStage=="Adult", SampleType=="skin swab") %>% nrow()
nrow(Kuen_toAdd_COMPLETE)

new_all_meta %>% filter(Host %in% multipleStudyHosts_COMPLETE) %>% filter(is.na(ENARUN), KUEN, LifeStage=="Adult", SampleType=="skin swab") %>%
  select(SampleID) %>% mutate(Full=TRUE) %>%
  full_join(Kuen_toAdd_COMPLETE %>% mutate(Full=FALSE)) 

# Kuen_toAdd_COMPLETE %>% View()
