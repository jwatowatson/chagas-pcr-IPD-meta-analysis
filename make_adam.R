library(tidyverse)
library(tibble)
library(dplyr)
library(kableExtra)
library(readr)

source('functions.R')

dm_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/DM 2024-04-15.csv')
ds_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/DS 2024-04-15.csv')
vs_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/VS 2024-04-15.csv')

sa_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/SA 2024-04-15.csv')

## need to automate unit checks
lb_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/LB 2024-04-15.csv')

in_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/IN 2024-04-15.csv')
ts_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/TS 2024-04-15.csv')
mb_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/MB 2024-04-15.csv')
pc_chagas = read_csv('CH_2024_1_Watson/DATA 2024-04-15/PC 2024-04-15.csv')

table(dm_chagas$STUDYID)


out=chagas_adam(dm_chagas = dm_chagas,
                ds_chagas = ds_chagas,
                in_chagas = in_chagas, 
                mb_chagas = mb_chagas,
                ts_chagas = ts_chagas,
                vs_chagas = vs_chagas,
                sa_chagas = sa_chagas,
                lb_chagas = lb_chagas,study_remove = 'CGBKZSR')


pcr_chagas=out$pcr_chagas
sa_chagas=out$sa_chagas
lb_chagas=out$lb_chagas

dat_summary = pcr_chagas %>% distinct(USUBJID, .keep_all = T)
table(dat_summary$ARM, dat_summary$BNZ_total_days)
# plot(pcr_chagas$VISIT_numeric, pcr_chagas$Day_frm_rand, col=as.numeric(as.factor(pcr_chagas$STUDYID)))
# legend('bottomright', col=1:3, legend = unique(pcr_chagas$STUDYID),lwd=2)



save(pcr_chagas, file = 'RData/pcr_chagas.RData')
save(sa_chagas, file = 'RData/sa_chagas.RData')
save(lb_chagas, file = 'RData/lb_chagas.RData')

# ### Try to get old ID variable
# 
# library(lubridate)
# load(file = 'RData/pcr_chagas.RData')
# pcr_chagas_all = pcr_chagas
# pcr_chagas = pcr_chagas %>% filter(VISIT=='Screening') %>% distinct(USUBJID,.keep_all = T)
# 
# 
# pcr_chagas$site_num = plyr::mapvalues(x = pcr_chagas$SITEID,
#                                       from = c('Cochabamba','Sucre','Tarija'),
#                                       to = c(1,3,2))
# table(pcr_chagas$site_num)
# pcr_chagas$year = year(as.Date(pcr_chagas$RFSTDTC,format = '%Y-%M'))
# 
# 
# dat = read_csv('~/Dropbox/MORU/ChagasEfficacy/Data/BENDITA/adser.csv')
# dat = dat %>% filter(!is.na(TRTP), VISIT=='Screening') %>% distinct(SUBJID,.keep_all = T)
# dat$ID = sapply(dat$SUBJID, function(x) unlist(strsplit(x,split = '-'))[2])
# unique(dat$TRTP)
# unique(pcr_chagas$ARM)
# 
# dat$year = year(dat$TDAT)
# 
# dat$ARM = 
#   plyr::mapvalues(x = dat$TRTP, 
#                   from = c("Placebo",
#                            "BZN 300 mg (Weekly) 8 wks / E1224 300 mg",
#                            "BZN 300 mg 4 wks",
#                            "BZN 150 mg 4 wks / E1224 300 mg",
#                            "BZN 300 mg 8 wks",
#                            "BZN 300 mg 2 wks",
#                            "BZN 150 mg 4 wks"),
#                   to = c('PLACEBO',                                     
#                          'BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE',
#                          'BENZNIDAZOLE 300MG 4W',                       
#                          'BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE',       
#                          'BENZNIDAZOLE 300MG 8W',                       
#                          'BENZNIDAZOLE 300MG 2W',                       
#                          'BENZNIDAZOLE 150MG 4W'))
# 
# 
# dat$SITEID = as.numeric(dat$SITEID)
# table(dat$SITEID)
# 
# dat$Unique_ID = apply(dat[, c('SEX','AGE','ARM','SITEID')],1,function(x) paste(x,collapse = '_'))
# pcr_chagas$Unique_ID = apply(pcr_chagas[, c('SEX','AGE','ARM','site_num')],1,function(x) paste(x,collapse = '_'))
# 
# all(sort(unique(dat$Unique_ID))==sort(unique(pcr_chagas$Unique_ID)))
# 
# 
# 
# 
# dat_orig1 = readxl::read_excel('CH_2024_1_Watson/DATA 2024-03-08/Planillas Revisadas - Proyecto Bendita (Oficial).xlsx',
#                                sheet = 'SLAN (S001-S143)')
# dat_orig2 = readxl::read_excel('CH_2024_1_Watson/DATA 2024-03-08/Planillas Revisadas - Proyecto Bendita (Oficial).xlsx',
#                                sheet = 'BIORAD (B001-B277)')
# 
# dat_orig_all = rbind(dat_orig1, dat_orig2) %>% filter(Fluor %in% c('Tc','FAM')) %>% 
#   ungroup() %>%
#   mutate(CT = ifelse(CT=='No Ct',40, CT),
#          CT = as.numeric(CT)) %>%
#   group_by(Label) %>%
#   mutate(ID = unlist(strsplit(Label[1],'-'))[1],
#          ID = as.numeric(ID),
#          VISIT = unlist(strsplit(Label[1],'-'))[2],
#          VISIT = ifelse(VISIT=='SCT','SCR',VISIT)) %>% 
#   filter(!is.na(ID), !is.na(VISIT), !is.na(CT), VISIT=='SCR') %>%
#   group_by(ID) %>%
#   mutate(
#     N_PCRS = length(CT),
#     N_pos = sum(CT<40)
#   ) %>% filter(ID %in% dat$ID) 
# 
# dat_orig = dat_orig_all %>% distinct(ID, .keep_all = T)
# 
# dat = merge(dat, dat_orig[, c('ID', 'N_pos')], by = 'ID', all = T)
# 
# pcr_chagas$ID_true = NA;k=1
# for(i in 1:nrow(pcr_chagas)){
#   id_unique = pcr_chagas$Unique_ID[i]
#   ind_match = which(dat$Unique_ID==id_unique)
#   if(length(ind_match)==1){
#     pcr_chagas$ID_true[i]=dat$ID[ind_match]
#     print(id_unique)
#     print(k)
#     k=k+1
#     dat = dat[-ind_match, ]
#   }
#   
# }
# 
# for(i in which(is.na(pcr_chagas$ID_true))){
#   id_unique = pcr_chagas$Unique_ID[i]
#   ind_match = which(
#     dat$Unique_ID==id_unique & 
#       !is.na(dat$N_pos) & 
#       dat$N_pos==pcr_chagas$N_pos[i])
#   
#     if(length(ind_match)==1){
#       pcr_chagas$ID_true[i]=dat$ID[ind_match]
#       dat = dat[-ind_match, ]
#     }
# }
# 
# View(pcr_chagas %>% filter(is.na(ID_true)))
# 
# pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='M_31_BENZNIDAZOLE 300MG 2W_2']='2062'
# pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_30_BENZNIDAZOLE 300MG 2W_3']='3088'
# pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_27_BENZNIDAZOLE 150MG 4W_3']='3031'
# pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_26_BENZNIDAZOLE 300MG 2W_2']='2019'
# pcr_chagas$ID_true[is.na(pcr_chagas$ID_true) & pcr_chagas$Unique_ID=='F_37_BENZNIDAZOLE 150MG 4W_2']='2053'
# 
# View(pcr_chagas %>% filter(is.na(ID_true)))
# 
# ID_map = pcr_chagas[, c('USUBJID','ID_true')]
# save(ID_map, file = 'RData/ID_map.RData')
# 
