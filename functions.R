# This function makes an ADAM dataset
chagas_adam = function(dm_chagas, ds_chagas, in_chagas, mb_chagas, ts_chagas,vs_chagas,sa_chagas, lb_chagas){
  
  # only select individuals who were randomised
  dm_chagas = dm_chagas %>% 
    filter(ARM != "NOT PROVIDED IN THE CONTRIBUTED DATASET") %>%
    mutate(
      ARM = gsub(pattern='BENZNIDAZOLE',replacement='BNZ',x=ARM),
      ARM = gsub(pattern='FOSRAVUCONAZOLE',replacement='E1224',x=ARM),
      ARM = case_when(
        ARM=='BNZ' ~ 'BNZ 300MG 8W',
        ARM=='FEXINIDAZOLE 3.6g' ~ 'FEX 1200MG 3D',
        ARM=='FEXINIDAZOLE 6.0g' ~ 'FEX 600MG 10D',
        ARM=='FEXINIDAZOLE 6.6g' ~ 'FEX 600MG 3D + 1200MG 4D',
        ARM=="E1224 SHORT DOSE 4w" ~ "E1224 400MG 4W",
        ARM=="E1224 HIGH DOSE 8w" ~ "E1224 400MG 8W",
        ARM=="E1224 LOW DOSE 8w" ~ "E1224 200MG 8W",
        ARM=="BNZ 300MG WEEKLY 8W E1224" ~ "BNZ 300MG/WK + E1224 300MG 8W",
        ARM=="BNZ 150MG 4W E1224" ~ "BNZ 150MG 4W + E1224 300MG 8W",
        T ~ ARM
      )
    )
  ds_chagas = ds_chagas %>% filter(DSTERM=='RANDOMIZED') %>%# extract the day of randomisation relative to screening
    mutate(Day_randomisation = DSSTDY)
  
  vs_chagas = vs_chagas%>%
    filter(USUBJID %in% dm_chagas$USUBJID,
           VSTEST=='Weight', 
           EPOCH%in%c('SCREENING','BASELINE'))%>%
    arrange(USUBJID)%>%group_by(USUBJID)%>%
    mutate(weight=mean(VSORRES))%>%distinct(USUBJID, .keep_all = T)%>%
    select(USUBJID, weight)
  
  dm_chagas = merge(dm_chagas, vs_chagas, by = 'USUBJID')
  
  # These two individuals have screening visits with an MBDY that is after their day 1 visit
  mb_chagas$MBDY[mb_chagas$USUBJID == '596' & mb_chagas$VISIT == 'Screening'] = 0
  mb_chagas$MBDY[mb_chagas$USUBJID == '603' & mb_chagas$VISIT == 'Screening'] = 0
  # get day of randomisation from DS and then compute days relative to randomisation
  mb_chagas = 
    merge(mb_chagas, ds_chagas[, c('USUBJID','Day_randomisation')], all = T) %>%
    mutate(
      Day_frm_rand = MBDY-Day_randomisation
    )
  
  in_chagas = merge(in_chagas, dm_chagas, by = c('USUBJID', 'STUDYID'))
  in_chagas_summary = 
    merge(in_chagas, ds_chagas[, c('USUBJID','Day_randomisation')], all = T) %>%
    mutate(
      IN_Day_frm_rand = INDY-Day_randomisation
    ) %>%
    filter(USUBJID %in% dm_chagas$USUBJID, !is.na(IN_Day_frm_rand)) %>%
    group_by(USUBJID) %>%
    filter(INTRT %in% c('FOSRAVUCONAZOLE','BENZNIDAZOLE','FEXINIDAZOLE','PLACEBO')) %>%
    mutate( # make dummy benznidazole entries, so that placebo ar recorded as zero doses of all drugs
      INDOSE = ifelse(INTRT=='PLACEBO',0,INDOSE),
      INTRT = ifelse(INTRT=='PLACEBO','BENZNIDAZOLE',INTRT)
    ) %>%
    group_by(USUBJID, IN_Day_frm_rand) %>%
    mutate(
      bnz_daily_dose = sum(INDOSE[INTRT=='BENZNIDAZOLE']),
      fos_daily_dose = sum(INDOSE[INTRT=='FOSRAVUCONAZOLE']),
      fex_daily_dose = sum(INDOSE[INTRT=='FEXINIDAZOLE']),
    ) %>% group_by(USUBJID) %>%
    distinct(USUBJID, INDY, .keep_all = T) %>% group_by(USUBJID) %>%
    mutate(
      BNZ_total_dose = sum(bnz_daily_dose),
      FOS_total_dose = sum(fos_daily_dose),
      FEX_total_dose = sum(fex_daily_dose),
      BNZ_total_dose = ifelse(is.na(BNZ_total_dose), 0, BNZ_total_dose),
      FOS_total_dose = ifelse(is.na(FOS_total_dose), 0, FOS_total_dose),
      FEX_total_dose = ifelse(is.na(FEX_total_dose), 0, FEX_total_dose),
      BNZ_total_days = ifelse(BNZ_total_dose==0, 0, max(IN_Day_frm_rand[which(bnz_daily_dose>0)])),
      FOS_total_days = ifelse(FOS_total_dose==0, 0, max(IN_Day_frm_rand[which(fos_daily_dose>0)])),
      FEX_total_days = ifelse(FEX_total_dose==0, 0, max(IN_Day_frm_rand[which(fex_daily_dose>0)]))
    ) %>%
    distinct(USUBJID, .keep_all = T) %>% 
    select(USUBJID, 
           starts_with('BNZ',ignore.case = F),
           starts_with('FOS',ignore.case = F),
           starts_with('FEX',ignore.case = F))
  
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
      CT = ifelse(CT>40, 40, CT)
    ) 
  
  
  pcr_chagas = merge(pcr_chagas, in_chagas_summary, by = 'USUBJID', all = T)
  pcr_chagas = merge(pcr_chagas, dm_chagas, by = c('USUBJID','STUDYID')) %>%
    mutate(
      STUDYID = case_when(
        STUDYID=='CGBKZSR' ~ 'FEX12',
        STUDYID=='CGLETTP' ~ 'E1224',
        STUDYID=='CGTNWOV' ~ 'BENDITA'),
      ##****not yet curated****##
      BNZ_total_dose = ifelse(STUDYID=='E1224' & ARM=='BNZ 300MG 8W', 16800, BNZ_total_dose),
      BNZ_total_days = ifelse(STUDYID=='E1224' & ARM=='BNZ 300MG 8W', 56, BNZ_total_days),
      BNZ_total_dose_mg_kg = BNZ_total_dose/weight,
      FOS_total_dose_mg_kg = FOS_total_dose/weight,
      FEX_total_dose_mg_kg = FEX_total_dose/weight
    )
  
  
  pcr_chagas$ID = apply(pcr_chagas[, c('STUDYID','USUBJID')], 1,
                        function(x) paste(c(x[1],as.numeric(x[2])), collapse ='_'))
  
  # Make end of treatment variable
  pcr_chagas = pcr_chagas %>%
    group_by(ID) %>%
    mutate(EOT = case_when(
      ARM == 'PLACEBO' ~ 0,
      ARM == 'BNZ 150MG 4W + E1224 300MG 8W' ~ 56,
      ARM == 'BNZ 150MG 4W' ~ 28,
      ARM == 'BNZ 300MG 4W' ~ 28,
      ARM == 'BNZ 300MG 8W' ~ 56,
      ARM == 'BNZ 300MG/WK + E1224 300MG 8W' ~ 56,
      ARM == 'BNZ 300MG 2W' ~ 14,
      ARM == 'FEX 1200MG 3D' ~ 3,
      ARM == 'FEX 600MG 10D' ~ 10,
      ARM == 'FEX 600MG 3D + 1200MG 4D' ~ 7,
      ARM == "E1224 400MG 8W" ~ 56,
      ARM == "E1224 200MG 8W" ~ 56,
      ARM == "E1224 400MG 4W" ~ 28
    )) %>%
    ungroup() %>%
    mutate(
      VISIT_numeric = as.numeric(gsub(pattern ='D',replacement = '',
                                      x = ifelse(VISIT_trans=='Screening',-1,VISIT_trans))))
  
  stopping_drug = sa_chagas %>% filter(SAACN %in% c("DRUG WITHDRAWN")) %>%
    group_by(USUBJID) %>%
    mutate(time_stop = min(SASTDY)) %>%
    distinct(USUBJID,.keep_all = T) %>% select(USUBJID, time_stop)
  
  pause_drug = sa_chagas %>% filter(SAACN %in% c("DRUG INTERRUPTED")) %>%
    group_by(USUBJID) %>%
    mutate(time_pause = min(SASTDY)) %>%
    distinct(USUBJID,.keep_all = T) %>% select(USUBJID, time_pause)
  
  pcr_chagas = merge(pcr_chagas, stopping_drug, by='USUBJID',all=T)
  pcr_chagas = merge(pcr_chagas, pause_drug, by='USUBJID',all=T)
  
  sa_chagas$SATERM = stringr::str_conv(sa_chagas$SATERM, "UTF-8")
  sa_chagas$SATERM = tolower(sa_chagas$SATERM)
  sa_chagas = sa_chagas %>% filter(EPOCH %in% c('FOLLOW-UP','TREATMENT') | 
                                     SACAT=='CONTRIBUTOR-REPORTED ADVERSE EVENTS',
                                   is.na(SAREL) | 
                                     SAREL %in% c('RELATED','UNASSESSABLE/UNCLASSIFIED'))
  sa_chagas =
    merge(sa_chagas, 
          dm_chagas, by = c('USUBJID','STUDYID')
    ) %>%
    mutate(
      STUDYID = case_when(
        STUDYID=='CGBKZSR' ~ 'FEX12',
        STUDYID=='CGLETTP' ~ 'E1224',
        STUDYID=='CGTNWOV' ~ 'BENDITA'))
  
  
  lb_chagas = lb_chagas %>% 
    mutate(
      LBORRES=as.numeric(LBORRES),
      LBORRES = ifelse(LBORRESU=='10^9/L', LBORRES*10^3, LBORRES),
      LBORRES = ifelse(LBORRESU=='g/L', LBORRES/10, LBORRES))
  
  lb_chagas_wide =
    merge(lb_chagas, 
          ds_chagas[, c('USUBJID','Day_randomisation')], 
          all = T,by = 'USUBJID') %>%
    mutate(
      Day_frm_rand = LBDY-Day_randomisation
    ) %>% 
    filter(!is.na(LBORRES), USUBJID %in% dm_chagas$USUBJID) %>% 
    pivot_wider(names_from = LBTEST, 
                values_from = LBORRES,
                id_cols = c('STUDYID','USUBJID','Day_frm_rand','VISIT'), 
                values_fn = mean)%>%
    mutate(
      STUDYID = case_when(
        STUDYID=='CGBKZSR' ~ 'FEX12',
        STUDYID=='CGLETTP' ~ 'E1224',
        STUDYID=='CGTNWOV' ~ 'BENDITA')) %>%
    arrange(STUDYID, USUBJID, Day_frm_rand)
  lb_chagas_wide =
    merge(lb_chagas_wide, 
          dm_chagas[, c('USUBJID','ARM')], 
          all = T,by = 'USUBJID')
  
  return(list(pcr_chagas=pcr_chagas, sa_chagas=sa_chagas, lb_chagas=lb_chagas_wide))
}



make_matrix_pcr = function(pcr_dat, N_max_PCR=9){
  
  require(RColorBrewer)
  if(is.na(N_max_PCR)) N_max_PCR=max(pcr_dat$PCR_number)
  pcr_dat = pcr_dat %>% filter(PCR_number<=N_max_PCR)
  
  K=length(unique(pcr_dat$VISIT_trans)) * N_max_PCR
  N=length(unique(pcr_dat$ID))
  
  pcr_dat$PCR_name = apply(pcr_dat[, c('PCR_number','VISIT_trans')], 1, function(x) paste(x[1],x[2], sep = '_'))
  pcr_mat = array(NA, dim=c(N,K))
  pcr_name_vals = pcr_dat %>% arrange(VISIT_numeric)
  
  colnames(pcr_mat) = apply(expand.grid(1:N_max_PCR, unique(pcr_name_vals$VISIT_trans)), 1, 
                            function(x) paste(x[1], x[2],sep = '_'))
  # print(colnames(pcr_mat))
  rownames(pcr_mat) = unique(pcr_dat$ID)
  
  for(id in unique(pcr_dat$ID)){
    ind = which(pcr_dat$ID==id)
    pcr_mat[id, pcr_dat$PCR_name[ind]] = pcr_dat$CT[ind]
  }
  return(pcr_mat)
  
}


plot_pcr_matrix = function(pcr_mat, my_breaks, my_break_legend, my_cols, 
                           visit_labels=NA, arm_labels, h_lines_ind, y_lab_ind,
                           cex_y_lab=1,
                           plot_legend=T){  
  
  arm_labels=gsub(pattern = '+',replacement = '\n+',x = arm_labels,fixed = T)
  
  image(t(pcr_mat), x=1:ncol(pcr_mat), y = 1:nrow(pcr_mat),
        breaks = my_breaks,
        col = my_cols,
        xlab='', xaxt = 'n', yaxt='n', ylab='')
  if(plot_legend) {
    legend('topright', inset= c(-0.12,0), title = 'CT',
           fill= my_cols, xpd=T,
           legend = my_break_legend)
  }
  if(is.na(visit_labels)){
    visit_labels = unique(unlist(sapply(colnames(pcr_mat), function(x) unlist(strsplit(x,split = '_'))[2])))
  }
  axis(1, at = (1:(ncol(pcr_mat)/9))*9 - 4.5, 
       labels = visit_labels, las=1,tick = F)
  
  abline(h=h_lines_ind, lty=2)
  
  abline(v=(1:(ncol(pcr_mat)/9))*9 + 0.5, lty=2)
  
  axis(2, at = y_lab_ind, labels = arm_labels,las=1,tick = F,cex.axis=cex_y_lab)
  
}



make_stan_dataset_model_v2 = function(pcr_chagas_ss,x_covs = NULL, PCR_spec=1,
                                      p_1_beta_prior=.5, p_2_beta_prior=2){
  
  pcr_chagas_ss = pcr_chagas_ss %>% arrange(ARM, ID, Day_frm_rand)
  pcr_mat = make_matrix_pcr(pcr_dat = pcr_chagas_ss,N_max_PCR = 9)
  
  par(mar=c(2,7,0.5,0.5))
  xx2 = pcr_chagas_ss %>% ungroup %>% distinct(ID,.keep_all = T)
  ind_y = which(!duplicated(xx2$ARM))
  plot_pcr_matrix(pcr_mat = pcr_mat,
                  my_breaks = my_breaks, 
                  my_break_legend = my_breaks_legend,
                  my_cols = my_cols,
                  arm_labels = unique(pcr_chagas_ss$ARM),
                  h_lines_ind = ind_y[-1]-0.5,
                  y_lab_ind = ind_y + diff(c(ind_y, nrow(pcr_mat)+1))/2,
                  cex_y_lab = .7,
                  plot_legend = F)
  
  pcr_summary = pcr_chagas_ss %>%
    distinct(ID, Day_frm_rand, MBREFID, .keep_all = T) %>%
    ungroup() %>% 
    mutate(
      ID_numeric = as.numeric(as.factor(ID)),
      ARM_numeric = as.numeric(as.factor(ARM))) %>%
    arrange(Mean_CT_inv_baseline, ID) %>% group_by(ID) %>%
    mutate(
      N_fup_samples = length(N_PCRs_sample)
    )
  pcr_unique = pcr_summary %>% distinct(ID, .keep_all = T)
  if(!is.null(x_covs)) {
    X_matrix = model.matrix(~. -1, data = pcr_unique[, x_covs, drop=F])
  } else {
    X_matrix = array(0, dim = c(nrow(pcr_unique),1))
  }
  table(pcr_unique$N_fup_samples, pcr_unique$ARM)
  ind_not_dup = which(!duplicated(pcr_summary$ID))
  data_list_stan = list(n_id = length(unique(pcr_summary$ID)), 
                        n_samples = nrow(pcr_summary),
                        Kmax = max(pcr_summary$N_PCRs_sample),
                        K_arms = length(unique(pcr_summary$ARM)),
                        trt_group = pcr_summary$ARM_numeric[ind_not_dup],
                        K_cov=ncol(X_matrix),
                        X=X_matrix,
                        y_pos = pcr_summary$N_pos_sample,
                        K_FUP = pcr_summary$N_PCRs_sample,
                        ind_start = ind_not_dup,
                        ind_end = c((ind_not_dup-1)[-1], nrow(pcr_summary)),
                        PCR_specificity=PCR_spec,
                        p_1_beta_prior=p_1_beta_prior,
                        p_2_beta_prior=p_2_beta_prior)
  
  return(data_list_stan)
}