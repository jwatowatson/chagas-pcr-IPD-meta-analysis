library(iddoverse)
library(tibble)
library(dplyr)
library(kableExtra)
library(readr)

dm_chagas = read_csv('CH_2024_1_Watson/DATA 2024-03-08/DM 2024-03-08.csv')
in_chagas = read_csv('CH_2024_1_Watson/DATA 2024-03-08/IN 2024-03-08.csv')

# select only patients who passed screening and were randomised
dm_chagas = dm_chagas %>% filter(ARM != "NOT PROVIDED IN THE CONTRIBUTED DATASET")
# sum dose by day
in_chagas = in_chagas %>% 
  filter(USUBJID %in% dm_chagas$USUBJID) %>%
  group_by(USUBJID) %>%
  mutate(
    Day1 = unique(INDY[!is.na(VISIT) & VISIT=='Day 1 (Post-Dose)']),
    Day_frm_rand = INDY - Day1) %>%
  filter(INTRT %in% c('FOSRAVUCONAZOLE','BENZNIDAZOLE')) %>%
  group_by(USUBJID, Day_frm_rand) %>%
  mutate(
    BNZ_daily_dose = sum(INTRT=='BENZNIDAZOLE')*150,
    FOS_daily_dose = sum(INTRT=='FOSRAVUCONAZOLE')*300
  ) %>% group_by(USUBJID) %>%
  distinct(USUBJID, Day_frm_rand, .keep_all = T)%>%
  mutate(
    BNZ_total_dose = sum(BNZ_daily_dose),
    FOS_total_dose = sum(FOS_daily_dose),
    BNZ_total_days = sum(BNZ_daily_dose>0),
    FOS_total_days = sum(FOS_daily_dose>0),
    Day_last_bnz = max(Day_frm_rand[BNZ_daily_dose>0]),
    Day_last_fos = ifelse(FOS_total_dose==0,NA,max(Day_frm_rand[FOS_daily_dose>0]))) %>%
  distinct(USUBJID, .keep_all = T)

ts_chagas = read_csv('CH_2024_1_Watson/DATA 2024-03-08/TS.csv')
mb_chagas = read_csv('CH_2024_1_Watson/DATA 2024-03-08/MB 2024-03-08.csv')

pcr_chagas = mb_chagas %>% filter(
  USUBJID %in% dm_chagas$USUBJID,
  MBMETHOD == "REAL-TIME POLYMERASE CHAIN REACTION ASSAY",
  MBTSTDTL == "QUANTIFICATION CYCLE NUMBER") %>% 
  mutate(MBORRES = ifelse(MBORRES=='>=40', 40, MBORRES),
         CT = as.numeric(MBORRES)) %>%
  group_by(USUBJID, VISIT) %>%
  mutate(Any_Pos_40 = any(CT<40),
         PCR_number = 1:length(CT),
         N_PCRs = max(PCR_number),
         N_pos = sum(CT<40)) %>%
  group_by(USUBJID, VISIT, MBGRPID) %>%
  mutate(N_PCRs_5ml_sample = length(CT),
         N_pos_5ml_sample = sum(CT<40))
pcr_chagas = merge(pcr_chagas, dm_chagas, by = c('USUBJID','STUDYID'))
pcr_chagas$ID = apply(pcr_chagas[, c('STUDYID','USUBJID')],1,function(x) paste(c(x[1],as.numeric(x[2])),collapse ='_'))

pcr_chagas = 
  merge(pcr_chagas, in_chagas[, c('USUBJID', "BNZ_daily_dose","Day_last_fos",
                                  "FOS_daily_dose","BNZ_total_dose","FOS_total_dose",
                                  "BNZ_total_days","FOS_total_days","Day_last_bnz")],
        by = 'USUBJID',all = T)

# Make relative day from randomisation
pcr_chagas$MBDY[pcr_chagas$USUBJID=='103'&pcr_chagas$VISIT=='Screening']=0
pcr_chagas$MBDY[pcr_chagas$USUBJID=='483'&pcr_chagas$VISIT=='Screening']=0

pcr_chagas = pcr_chagas %>% 
  group_by(ID) %>%
  mutate(
    ref_day = ifelse(sum(VISIT=='Day 2')>0, 'Day 2', 'Day 3'),
    day_0 = ifelse(ref_day=='Day 2', MBDY[VISIT=='Day 2']-2, MBDY[VISIT=='Day 3']-3),
    Day_frm_rand = MBDY-day_0)
# View(pcr_chagas %>% filter(Day_frm_rand==58, VISIT=='Week 10'))


pcr_chagas=pcr_chagas %>%
  mutate(ARM = factor(ARM, 
                      levels=c("PLACEBO",
                               "BENZNIDAZOLE 300MG 2W",
                               "BENZNIDAZOLE 150MG 4W", 
                               "BENZNIDAZOLE 300MG 4W",
                               "BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE",
                               "BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE",
                               "BENZNIDAZOLE 300MG 8W"))) %>%
  # filter(ID != 'CGTNWOV_347') %>% #person who only took 3 days
  group_by(ID) %>%
  mutate(EOT = case_when(
    ARM=='PLACEBO' ~ 0,
    ARM=='BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE' ~ 56,
    ARM=='BENZNIDAZOLE 150MG 4W' ~ 28,
    ARM=='BENZNIDAZOLE 300MG 4W' ~ 28,
    ARM=='BENZNIDAZOLE 300MG 8W' ~ 56,
    ARM=='BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE' ~ 56,
    ARM=='BENZNIDAZOLE 300MG 2W' ~ 14
  ))



save(pcr_chagas, file = 'RData/pcr_chagas.RData')
### Try to get old ID variable

library(lubridate)
load(file = 'RData/pcr_chagas.RData')
pcr_chagas_all = pcr_chagas
pcr_chagas = pcr_chagas %>% filter(VISIT=='Screening') %>% distinct(USUBJID,.keep_all = T)


pcr_chagas$site_num = plyr::mapvalues(x = pcr_chagas$SITEID,
                                      from = c('Cochabamba','Sucre','Tarija'),
                                      to = c(1,3,2))
table(pcr_chagas$site_num)
pcr_chagas$year = year(as.Date(pcr_chagas$RFSTDTC,format = '%Y-%M'))


dat = read_csv('~/Dropbox/MORU/ChagasEfficacy/Data/BENDITA/adser.csv')
dat = dat %>% filter(!is.na(TRTP), VISIT=='Screening') %>% distinct(SUBJID,.keep_all = T)
dat$ID = sapply(dat$SUBJID, function(x) unlist(strsplit(x,split = '-'))[2])
unique(dat$TRTP)
unique(pcr_chagas$ARM)

dat$year = year(dat$TDAT)

dat$ARM = 
  plyr::mapvalues(x = dat$TRTP, 
                  from = c("Placebo",
                           "BZN 300 mg (Weekly) 8 wks / E1224 300 mg",
                           "BZN 300 mg 4 wks",
                           "BZN 150 mg 4 wks / E1224 300 mg",
                           "BZN 300 mg 8 wks",
                           "BZN 300 mg 2 wks",
                           "BZN 150 mg 4 wks"),
                  to = c('PLACEBO',                                     
                         'BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE',
                         'BENZNIDAZOLE 300MG 4W',                       
                         'BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE',       
                         'BENZNIDAZOLE 300MG 8W',                       
                         'BENZNIDAZOLE 300MG 2W',                       
                         'BENZNIDAZOLE 150MG 4W'))


dat$SITEID = as.numeric(dat$SITEID)
table(dat$SITEID)

dat$Unique_ID = apply(dat[, c('SEX','AGE','ARM','SITEID')],1,function(x) paste(x,collapse = '_'))
pcr_chagas$Unique_ID = apply(pcr_chagas[, c('SEX','AGE','ARM','site_num')],1,function(x) paste(x,collapse = '_'))

all(sort(unique(dat$Unique_ID))==sort(unique(pcr_chagas$Unique_ID)))




dat_orig1 = readxl::read_excel('CH_2024_1_Watson/DATA 2024-03-08/Planillas Revisadas - Proyecto Bendita (Oficial).xlsx',
                               sheet = 'SLAN (S001-S143)')
dat_orig2 = readxl::read_excel('CH_2024_1_Watson/DATA 2024-03-08/Planillas Revisadas - Proyecto Bendita (Oficial).xlsx',
                               sheet = 'BIORAD (B001-B277)')

dat_orig_all = rbind(dat_orig1, dat_orig2) %>% filter(Fluor %in% c('Tc','FAM')) %>% 
  ungroup() %>%
  mutate(CT = ifelse(CT=='No Ct',40, CT),
         CT = as.numeric(CT)) %>%
  group_by(Label) %>%
  mutate(ID = unlist(strsplit(Label[1],'-'))[1],
         ID = as.numeric(ID),
         VISIT = unlist(strsplit(Label[1],'-'))[2],
         VISIT = ifelse(VISIT=='SCT','SCR',VISIT)) %>% 
  filter(!is.na(ID), !is.na(VISIT), !is.na(CT), VISIT=='SCR') %>%
  group_by(ID) %>%
  mutate(
    N_PCRS = length(CT),
    N_pos = sum(CT<40)
  ) %>% filter(ID %in% dat$ID) 

dat_orig = dat_orig_all %>% distinct(ID, .keep_all = T)

dat = merge(dat, dat_orig[, c('ID', 'N_pos')], by = 'ID', all = T)

pcr_chagas$ID_true = NA;k=1
for(i in 1:nrow(pcr_chagas)){
  id_unique = pcr_chagas$Unique_ID[i]
  ind_match = which(dat$Unique_ID==id_unique)
  if(length(ind_match)==1){
    pcr_chagas$ID_true[i]=dat$ID[ind_match]
    print(id_unique)
    print(k)
    k=k+1
    dat = dat[-ind_match, ]
  }
  
}

for(i in which(is.na(pcr_chagas$ID_true))){
  id_unique = pcr_chagas$Unique_ID[i]
  ind_match = which(
    dat$Unique_ID==id_unique & 
      !is.na(dat$N_pos) & 
      dat$N_pos==pcr_chagas$N_pos[i])
  
    if(length(ind_match)==1){
      pcr_chagas$ID_true[i]=dat$ID[ind_match]
      dat = dat[-ind_match, ]
    }
}

View(pcr_chagas %>% filter(is.na(ID_true)))

pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='M_31_BENZNIDAZOLE 300MG 2W_2']='2062'
pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_30_BENZNIDAZOLE 300MG 2W_3']='3088'
pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_27_BENZNIDAZOLE 150MG 4W_3']='3031'
pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_26_BENZNIDAZOLE 300MG 2W_2']='2019'
pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_37_BENZNIDAZOLE 150MG 4W_2']='2053'

View(pcr_chagas %>% filter(is.na(ID_true)))

ID_map = pcr_chagas[, c('USUBJID','ID_true')]
save(ID_map, file = 'RData/ID_map.RData')

