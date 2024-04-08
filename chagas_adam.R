chagas_adam = function(dm_chagas, in_chagas, mb_chagas, ts_chagas){
  dm_chagas = dm_chagas %>% 
    filter(ARM != "NOT PROVIDED IN THE CONTRIBUTED DATASET")
  
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
  
  pcr_chagas = mb_chagas %>%
    filter(
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
  pcr_chagas$ID = apply(pcr_chagas[, c('STUDYID','USUBJID')], 1,
                        function(x) paste(c(x[1],as.numeric(x[2])), collapse ='_'))
  
  pcr_chagas = 
    merge(pcr_chagas, in_chagas[, c('USUBJID', "BNZ_daily_dose","Day_last_fos",
                                    "FOS_daily_dose","BNZ_total_dose","FOS_total_dose",
                                    "BNZ_total_days","FOS_total_days","Day_last_bnz")],
          by = 'USUBJID', all = T)
  
  # Make relative day from randomisation
  pcr_chagas$MBDY[pcr_chagas$USUBJID == '103' & pcr_chagas$VISIT == 'Screening'] = 0
  pcr_chagas$MBDY[pcr_chagas$USUBJID == '483' & pcr_chagas$VISIT == 'Screening'] = 0
  
  pcr_chagas = pcr_chagas %>% 
    group_by(ID) %>%
    mutate(
      ref_day = ifelse(sum(VISIT == 'Day 2') > 0, 'Day 2', 'Day 3'),
      day_0 = ifelse(ref_day == 'Day 2', MBDY[VISIT == 'Day 2'] - 2, MBDY[VISIT == 'Day 3'] - 3),
      Day_frm_rand = MBDY-day_0)
  
  
  pcr_chagas = pcr_chagas %>%
    mutate(ARM = factor(ARM, 
                        levels = c("PLACEBO",
                                 "BENZNIDAZOLE 300MG 2W",
                                 "BENZNIDAZOLE 150MG 4W", 
                                 "BENZNIDAZOLE 300MG 4W",
                                 "BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE",
                                 "BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE",
                                 "BENZNIDAZOLE 300MG 8W"))) %>%
    # filter(ID != 'CGTNWOV_347') %>% #person who only took 3 days
    group_by(ID) %>%
    mutate(EOT = case_when(
      ARM == 'PLACEBO' ~ 0,
      ARM == 'BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE' ~ 56,
      ARM == 'BENZNIDAZOLE 150MG 4W' ~ 28,
      ARM == 'BENZNIDAZOLE 300MG 4W' ~ 28,
      ARM == 'BENZNIDAZOLE 300MG 8W' ~ 56,
      ARM == 'BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE' ~ 56,
      ARM == 'BENZNIDAZOLE 300MG 2W' ~ 14
    ))
  
  return(pcr_chagas)
}
