# This function makes an ADAM dataset
chagas_adam = function(dm_chagas, in_chagas, mb_chagas, ts_chagas){
  
  # only select individuals who were randomised
  dm_chagas = dm_chagas %>% 
    filter(ARM != "NOT PROVIDED IN THE CONTRIBUTED DATASET") %>%
    mutate(
      ARM = case_when(
        ARM=='BENZNIDAZOLE' ~ 'BENZNIDAZOLE 300MG 8W',
        ARM=='FEXINIDAZOLE 3.6g' ~ 'FEX 1200MG 3D',
        ARM=='FEXINIDAZOLE 6.0g' ~ 'FEX 600MG 10D',
        ARM=='FEXINIDAZOLE 6.6g' ~ 'FEX 600MG 3D + 1200MG 4D',
        T ~ ARM
      )
    )
  
  ## The solution for getting days from randomisation is very ad hoc and not great...
  in_chagas = in_chagas %>% 
    filter(USUBJID %in% dm_chagas$USUBJID, !is.na(INDY)) %>%
    group_by(USUBJID) %>%
    filter(INTRT %in% c('FOSRAVUCONAZOLE','BENZNIDAZOLE','FEXINIDAZOLE','PLACEBO')) %>%
    group_by(USUBJID, INDY) %>%
    mutate(
      BNZ_daily_dose = sum(INTRT=='BENZNIDAZOLE')*150,
      FOS_daily_dose = sum(INTRT=='FOSRAVUCONAZOLE')*300
    ) %>% group_by(USUBJID) %>%
    distinct(USUBJID, INDY, .keep_all = T)%>%
    mutate(
      BNZ_total_dose = sum(BNZ_daily_dose),
      FOS_total_dose = sum(FOS_daily_dose),
      BNZ_total_days = sum(BNZ_daily_dose>0),
      FOS_total_days = sum(FOS_daily_dose>0),
      # Day_last_bnz = max(INDY[BNZ_daily_dose>0]),
      # Day_last_fos = ifelse(FOS_total_dose==0,NA,max(INDY[FOS_daily_dose>0])
      ) %>%
    distinct(USUBJID, .keep_all = T)
  
  pcr_chagas = mb_chagas %>%
    filter(
      USUBJID %in% dm_chagas$USUBJID,
      MBMETHOD == "REAL-TIME POLYMERASE CHAIN REACTION ASSAY",
      MBTSTDTL == "QUANTIFICATION CYCLE NUMBER") %>% 
    mutate(
      VISIT_trans = case_when( ## tried to make consistent as nominal days at start of week
        VISIT=='Screening' ~ 'Screening',
        VISIT=="Day 1 (Post-treatment)" ~ 'D1',
        VISIT=="Day 1 (Post-Dose)" ~ 'D1',
        VISIT=="Day 2" ~ 'D2',
        VISIT=="Day 3" ~ 'D3',
        VISIT=="Week 2 - Day 8" ~ 'D7',
        VISIT=="D8" ~ 'D7',
        VISIT=='Week 2' ~ 'D7',
        VISIT=="Week 3 - Day 15" ~ 'D14',
        VISIT=='Week 3' ~ 'D14',
        VISIT=="D15" ~ 'D14',
        VISIT=='Week 4 - Day 22' ~ 'D21',
        VISIT=='Week 4' ~ 'D21',
        VISIT=='Week 6 - Day 36' ~ 'D35',
        VISIT=='Week 6' ~ 'D35',
        VISIT=='D36' ~ 'D35',
        VISIT=='Week 10 - Day 64' ~ 'D63',
        VISIT=='Week 10' ~ 'D63',
        VISIT=='D65' ~ 'D63',
        VISIT=='Week 12' ~ 'D77',
        VISIT=='Week 12 - Day 90' ~ 'D77',
        VISIT=='4 Month (Week 18)' ~ 'D120',
        VISIT=='Month 4 - Day 120' ~ 'D120',
        VISIT=='4M' ~ 'D120',
        VISIT=='6 Month (Week 27)' ~ 'D180',
        VISIT=='6M' ~ 'D180',
        VISIT=='Month 6 - Day 183' ~ 'D180',
        VISIT=='12 Month (Week 54)' ~ 'D360',
        VISIT=='Month 12 - Day 360' ~ 'D360',
        VISIT=='12M' ~ 'D360',
        T ~ NA
      ),
      MBORRES = ifelse(MBORRES=='>=40', 40, MBORRES),
      CT = as.numeric(MBORRES),
      CT = ifelse(CT>40, 40, CT)) %>%
    group_by(USUBJID) %>%
    mutate(
      Ref_day = case_when(
        any(VISIT_trans == 'D1') ~ 'D1',
        any(VISIT_trans == 'D2') ~ 'D2',
        any(VISIT_trans == 'D7') ~ 'D7',
        any(VISIT_trans == 'D14') ~ 'D14',
        T ~ NA
      ))%>% filter(!is.na(Ref_day)) %>%
    mutate(
      Ref_MBDY = unique(MBDY[VISIT_trans==Ref_day]),
      Day_frm_rand = MBDY - Ref_MBDY,
      Day_frm_rand = 
        case_when(
          Ref_day=='D1' ~ Day_frm_rand,
          Ref_day=='D2' ~ Day_frm_rand+1,
          Ref_day=='D7' ~ Day_frm_rand+6,
          Ref_day=='D14' ~ Day_frm_rand+13
        )
    )%>%
    group_by(USUBJID, VISIT_trans) %>%
    mutate(Any_Pos_40 = any(CT<40),
           PCR_number = 1:length(CT),
           N_PCRs = max(PCR_number),
           N_pos = sum(CT<40)) %>%
    group_by(USUBJID, VISIT, MBREFID) %>%
    mutate(N_PCRs_5ml_sample = length(CT),
           N_pos_5ml_sample = sum(CT<40))
  
  pcr_chagas = merge(pcr_chagas, dm_chagas, by = c('USUBJID','STUDYID')) %>%
   mutate(STUDYID = case_when(
      STUDYID=='CGBKZSR' ~ 'FEX12',
      STUDYID=='CGLETTP' ~ 'E1224',
      STUDYID=='CGTNWOV' ~ 'BENDITA',
    ))
  
  pcr_chagas$ID = apply(pcr_chagas[, c('STUDYID','USUBJID')], 1,
                        function(x) paste(c(x[1],as.numeric(x[2])), collapse ='_'))
  
  # pcr_chagas = 
  #   merge(pcr_chagas, in_chagas[, c('USUBJID', "BNZ_daily_dose","Day_last_fos",
  #                                   "FOS_daily_dose","BNZ_total_dose","FOS_total_dose",
  #                                   "BNZ_total_days","FOS_total_days","Day_last_bnz")],
  #         by = 'USUBJID', all = T)
  
  # Make relative day from randomisation
  # These two individuals have screening visits with an MBDY that is after their day 1 visit
  pcr_chagas$MBDY[pcr_chagas$USUBJID == '596' & pcr_chagas$VISIT == 'Screening'] = 0
  pcr_chagas$MBDY[pcr_chagas$USUBJID == '603' & pcr_chagas$VISIT == 'Screening'] = 0
  
  pcr_chagas = pcr_chagas %>%
    # mutate(ARM = factor(ARM, 
    #                     levels = c("PLACEBO",
    #                                "BENZNIDAZOLE 300MG 2W",
    #                                "BENZNIDAZOLE 150MG 4W", 
    #                                "BENZNIDAZOLE 300MG 4W",
    #                                "BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE",
    #                                "BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE",
    #                                "BENZNIDAZOLE 300MG 8W",
    #                                "FEXINIDAZOLE 3.6g",
    #                                "FEXINIDAZOLE 6.0g",
    #                                "FEXINIDAZOLE 6.6g",
    #                                "FOSRAVUCONAZOLE HIGH DOSE 8w",
    #                                "FOSRAVUCONAZOLE LOW DOSE 8w",
    #                                "FOSRAVUCONAZOLE SHORT DOSE 4w"))) %>%
    # filter(ID != 'CGTNWOV_347') %>% #person who only took 3 days
    group_by(ID) %>%
    mutate(EOT = case_when(
      ARM == 'PLACEBO' ~ 0,
      ARM == 'BENZNIDAZOLE 150MG 4W FOSRAVUCONAZOLE' ~ 56,
      ARM == 'BENZNIDAZOLE 150MG 4W' ~ 28,
      ARM == 'BENZNIDAZOLE 300MG 4W' ~ 28,
      ARM == 'BENZNIDAZOLE 300MG 8W' ~ 56,
      ARM == 'BENZNIDAZOLE 300MG WEEKLY 8W FOSRAVUCONAZOLE' ~ 56,
      ARM == 'BENZNIDAZOLE 300MG 2W' ~ 14,
      ARM == 'FEX 1200MG 3D' ~ 3,
      ARM == 'FEX 600MG 10D' ~ 10,
      ARM == 'FEX 600MG 3D + 1200MG 4D' ~ 7,
      ARM == "FOSRAVUCONAZOLE HIGH DOSE 8w" ~ 56,
      ARM == "FOSRAVUCONAZOLE LOW DOSE 8w" ~ 56,
      ARM == "FOSRAVUCONAZOLE SHORT DOSE 4w" ~ 28
    ))
  
  return(pcr_chagas)
}



make_matrix_pcr = function(pcr_dat, N_max_PCR=NA){
  
  require(RColorBrewer)
  if(is.na(N_max_PCR)) N_max_PCR=max(pcr_dat$PCR_number)
  pcr_dat = pcr_dat %>% filter(PCR_number<=N_max_PCR)
  
  K=length(unique(pcr_dat$VISIT_trans)) * N_max_PCR
  N=length(unique(pcr_dat$ID))
  
  pcr_dat$PCR_name = apply(pcr_dat[, c('PCR_number','VISIT_trans')], 1, function(x) paste(x[1],x[2], sep = '_'))
  pcr_mat = array(NA, dim=c(N,K))
  pcr_name_vals = pcr_dat %>% arrange(Day_frm_rand)
  
  colnames(pcr_mat) = apply(expand.grid(1:N_max_PCR, unique(pcr_name_vals$VISIT_trans)), 1, 
                            function(x) paste(x[1], x[2],sep = '_'))
  print(colnames(pcr_mat))
  rownames(pcr_mat) = unique(pcr_dat$ID)
  
  for(id in unique(pcr_dat$ID)){
    ind = which(pcr_dat$ID==id)
    pcr_mat[id, pcr_dat$PCR_name[ind]] = pcr_dat$CT[ind]
  }
  return(pcr_mat)
  
}
plot_pcr_matrix = function(pcr_mat, my_breaks, my_break_legend, my_cols, 
                           visit_labels, arm_labels,
                           ind_y){  
  image(t(pcr_mat), x=1:ncol(pcr_mat), y = 1:nrow(pcr_mat),
        breaks = my_breaks,
        col = my_cols,
        xlab='', xaxt = 'n', yaxt='n', ylab='')
  legend('topright', inset= c(-0.1,0), title = 'CT',
         fill= my_cols, xpd=T,
         legend = my_break_legend)
  axis(1, at = (1:(ncol(pcr_mat)/9))*9 - 4.5, 
       labels = visit_labels, las=2)
  
  abline(h=ind_y[-1], lty=2)
  
  abline(v=(1:(ncol(pcr_mat)/9))*9 + 0.5, lty=2)
  
  axis(2, at = ind_y+15, labels = arm_labels,las=1,tick = F)
  
}