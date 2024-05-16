# This function makes an ADAM dataset
chagas_adam = function(dm_chagas, ds_chagas, 
                       in_chagas, mb_chagas, 
                       ts_chagas, vs_chagas,
                       sa_chagas, lb_chagas,
                       study_remove=NULL){
  
  # only select individuals who were randomised
  dm_chagas = dm_chagas %>% 
    filter(#ARM != "NOT PROVIDED IN THE CONTRIBUTED DATASET",
      !STUDYID %in% study_remove) %>%
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
        ARM=='NOT PROVIDED IN THE CONTRIBUTED DATASET' ~ 'FAILED SCREENING',
        T ~ ARM
      )
    )
  
  ds_chagas = ds_chagas %>% filter(DSTERM=='RANDOMIZED',!STUDYID %in% study_remove) %>%# extract the day of randomisation relative to screening
    mutate(Day_randomisation = DSSTDY)
  
  vs_chagas = vs_chagas%>%
    filter(USUBJID %in% dm_chagas$USUBJID,
           !STUDYID %in% study_remove,
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
  mb_chagas = mb_chagas %>% filter(!STUDYID %in% study_remove)
  mb_chagas = 
    merge(mb_chagas, ds_chagas[, c('USUBJID','Day_randomisation')], all = T) %>%
    mutate(
      Day_frm_rand = MBDY-Day_randomisation
    )
  
  in_chagas = in_chagas %>% filter(!STUDYID %in% study_remove)
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
      # fex_daily_dose = sum(INDOSE[INTRT=='FEXINIDAZOLE']),
    ) %>% group_by(USUBJID) %>%
    distinct(USUBJID, INDY, .keep_all = T) %>% group_by(USUBJID) %>%
    mutate(
      BNZ_total_dose = sum(bnz_daily_dose),
      FOS_total_dose = sum(fos_daily_dose),
      # FEX_total_dose = sum(fex_daily_dose),
      BNZ_total_dose = ifelse(is.na(BNZ_total_dose), 0, BNZ_total_dose),
      FOS_total_dose = ifelse(is.na(FOS_total_dose), 0, FOS_total_dose),
      # FEX_total_dose = ifelse(is.na(FEX_total_dose), 0, FEX_total_dose),
      BNZ_total_days = ifelse(BNZ_total_dose==0, 0, max(IN_Day_frm_rand[which(bnz_daily_dose>0)])),
      FOS_total_days = ifelse(FOS_total_dose==0, 0, max(IN_Day_frm_rand[which(fos_daily_dose>0)])),
      # FEX_total_days = ifelse(FEX_total_dose==0, 0, max(IN_Day_frm_rand[which(fex_daily_dose>0)]))
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
      # FEX_total_dose_mg_kg = FEX_total_dose/weight
    )
  
  # Make end of treatment variable
  pcr_chagas = pcr_chagas %>%
    group_by(USUBJID) %>%
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
  
  stopping_drug = sa_chagas %>% filter(SAACN %in% c("DRUG WITHDRAWN"), !STUDYID %in% study_remove) %>%
    group_by(USUBJID) %>%
    mutate(time_stop = min(SASTDY)) %>%
    distinct(USUBJID,.keep_all = T) %>% select(USUBJID, time_stop)
  
  pause_drug = sa_chagas %>% filter(SAACN %in% c("DRUG INTERRUPTED"), !STUDYID %in% study_remove) %>%
    group_by(USUBJID) %>%
    mutate(time_pause = min(SASTDY)) %>%
    distinct(USUBJID,.keep_all = T) %>% select(USUBJID, time_pause)
  
  pcr_chagas = merge(pcr_chagas, stopping_drug, by='USUBJID',all=T)
  pcr_chagas = merge(pcr_chagas, pause_drug, by='USUBJID',all=T)
  
  sa_chagas$SATERM = stringr::str_conv(sa_chagas$SATERM, "UTF-8")
  sa_chagas$SATERM = tolower(sa_chagas$SATERM)
  sa_chagas2 = sa_chagas %>% filter(EPOCH %in% c('FOLLOW-UP','TREATMENT') | 
                                     SACAT=='CONTRIBUTOR-REPORTED ADVERSE EVENTS',
                                   is.na(SAREL) | 
                                     SAREL %in% c('RELATED','UNASSESSABLE/UNCLASSIFIED'),
                                   !STUDYID %in% study_remove)
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
  N=length(unique(pcr_dat$USUBJID))
  
  ## order correctly
  pcr_dat = pcr_dat %>% 
    arrange(ARM, Percent_pos_baseline, USUBJID, Day_frm_rand, MBREFID, MBGRPID) %>%
    group_by(USUBJID, VISIT_trans) %>%
    mutate(PCR_number= as.character(1:n()))
  
  pcr_dat$PCR_name = apply(pcr_dat[, c('PCR_number','VISIT_trans')], 1, function(x) paste((x[1]),x[2], sep = '_'))
  pcr_mat = array(NA, dim=c(N,K))
  pcr_name_vals = pcr_dat %>% arrange(VISIT_numeric)
  
  colnames(pcr_mat) = apply(expand.grid(1:N_max_PCR, unique(pcr_name_vals$VISIT_trans)), 1, 
                            function(x) paste(x[1], x[2],sep = '_'))
  # print(colnames(pcr_mat))
  rownames(pcr_mat) = unique(pcr_dat$USUBJID)
  
  for(id in unique(pcr_dat$USUBJID)){
    ind = which(pcr_dat$USUBJID==id)
    pcr_mat[as.character(id), pcr_dat$PCR_name[ind]] = pcr_dat$CT[ind]
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



make_stan_dataset_model_v2 = 
  function(pcr_chagas_ss,x_covs = NULL, 
           PCR_spec=1,
           p_1_beta_prior=.5, 
           p_2_beta_prior=2){
    
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
    
    print(cbind(unique(pcr_summary$ARM), unique(pcr_summary$ARM_numeric)))
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




make_stan_dataset_model_quant = 
  function(pcr_chagas_ss_input,
           p_1_beta_prior=.5,lambda_max=10,
           plot_data=T)
  {
    
    pcr_chagas_ss = pcr_chagas_ss_input %>%
      group_by(ID_numeric, Day_frm_rand, MBREFID) %>%
      mutate(
        # unique Blood Sample ID
        Sample_ID = paste(c(unique(ID_numeric),unique(Day_frm_rand),unique(MBREFID)),collapse = '_'),
        CT_mean = mean(CT) # Mean CT value from that blood sample
      ) %>%
      group_by(ID_numeric, Day_frm_rand) %>%
      mutate(
        # make a unique Timepoint ID
        Timepoint_ID = paste(c(unique(ID_numeric),unique(Day_frm_rand)),collapse = '_')
      ) %>%
      select(ID_numeric,USUBJID,Sample_ID,MBREFID,Baseline_CT,VISIT_trans,VISIT_numeric,
             VISIT,CT,Day_frm_rand,ARM,EOT,PCR_number,CT_mean,Percent_pos_baseline,
             N_samples, N_PCRs_sample,Timepoint_ID,Sample_No,N_pos, N_PCRs, MBGRPID)
    
    xx2 = pcr_chagas_ss %>% ungroup %>% distinct(ID_numeric,.keep_all = T)
    trt = as.numeric(xx2$ARM)
    
    if(plot_data){
      pcr_mat = make_matrix_pcr(pcr_dat = pcr_chagas_ss,N_max_PCR = 9)
      # plot data used for model fitting
      par(mar=c(2,7,0.5,0.5))
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
    }
    # select unique timepoints
    pcr_timepoints = pcr_chagas_ss %>% 
      distinct(ID_numeric, Day_frm_rand,.keep_all = T)
    
    treatment = as.numeric(pcr_timepoints$Day_frm_rand>pcr_timepoints$EOT)
    print(table(pcr_timepoints$ARM, treatment))
    
    Smax=max(pcr_chagas_ss$N_samples)
    Kmax=max(pcr_chagas_ss$N_PCRs_sample)
    
    # number of blood samples per timepoint
    j_samples = pcr_timepoints$N_samples
    
    # number of PCR replicates per blood sample
    k_replicates=pcr_chagas_ss %>% 
      pivot_wider(id_cols = c(ID_numeric, Day_frm_rand), 
                  names_from = Sample_No, 
                  values_from = N_PCRs_sample, values_fn = max) %>%
      ungroup() %>%
      select(-ID_numeric, -Day_frm_rand)
    k_replicates[is.na(k_replicates)]=0
    k_replicates=as.matrix(k_replicates)
    
    id_ind = pcr_timepoints$ID_numeric
    ind_start = which(!duplicated(pcr_timepoints$ID_numeric))
    ind_end = as.integer(c((ind_start-1)[-1], nrow(pcr_timepoints)))
    
    CT_obs = array(40, dim = c(nrow(pcr_timepoints),Smax,Kmax))
    rownames(CT_obs) = pcr_timepoints$Timepoint_ID
    
    for(t in 1:nrow(pcr_timepoints)){
      for(j in 1:j_samples[t]){
        ind = which(pcr_chagas_ss$Timepoint_ID==rownames(CT_obs)[t] &
                      pcr_chagas_ss$Sample_No==j)
        if(length(ind)>0){
          CT_obs[t, j, 1:length(ind)] = pcr_chagas_ss$CT[ind]
        }
      }
    }
    
    # get summaries of number of + PCRs per timepoint
    y_pos = pcr_timepoints$N_pos
    K_total = pcr_timepoints$N_PCRs
    #check number of positives is less than or equal to the number of PCRs
    !any(y_pos > K_total)
    
    # record day from randomisation for each timepoint
    t_actual = pcr_timepoints$Day_frm_rand
    
    data_list_stan = list(N = length(trt), 
                          Tmax = nrow(pcr_timepoints),
                          Smax=Smax,
                          Kmax = Kmax,
                          trt = as.array(trt),
                          Trt_max=max(trt),
                          EOT = treatment,
                          j_samples=j_samples,
                          k_replicates=k_replicates,
                          CT_obs = CT_obs,
                          y_pos = y_pos,
                          K_total = K_total,
                          id_ind = as.array(id_ind),
                          ind_start = as.array(ind_start),
                          ind_end = as.array(ind_end),
                          p_1_beta_prior=p_1_beta_prior,
                          lambda_max = lambda_max,
                          baseline_CT = xx2$Baseline_CT,
                          t_actual=as.array(t_actual),
                          CT_blood_sample_mean = pcr_timepoints$CT_mean,
                          pcr_chagas_ss=pcr_chagas_ss)
    
    return(data_list_stan)
  }



get_analysis_dataset = function(dat, arms_select=NA, analysis_IDs=NA, 
                                max_PCR_replicates=3,
                                EOT_delta=0,
                                EOT_placebo = 0){
  
  # Select IDs, timepoints and make numeric ID for stan code
  if(any(is.na(arms_select))) arms_select=unique(dat$ARM)
  if(any(is.na(analysis_IDs))) analysis_IDs=unique(dat$USUBJID)
  
  pcr_chagas_ss_input = dat %>% ungroup() %>%
    filter(ARM %in% arms_select,
           USUBJID %in% analysis_IDs,
           Day_frm_rand <=0 | #baseline
             Day_frm_rand > EOT+EOT_delta | #follow-up
             ARM=='PLACEBO' # not treated
    ) %>%
    mutate(EOT = ifelse(ARM=='PLACEBO', EOT_placebo, EOT)) %>%
    arrange(ARM, USUBJID) %>%
    group_by(ARM,USUBJID) %>%
    mutate(
      ID_numeric = cur_group_id() # make a numeric ID 1..N
    ) %>%
    group_by(ID_numeric, Day_frm_rand) %>%
    mutate(
      N_samples = length(unique(MBREFID)), # no of blood samples
      Sample_No = as.numeric(as.factor(MBREFID)),
      PCR_number = 1:length(CT) # no of PCRs
    ) %>%
    group_by(ID_numeric, VISIT, MBREFID) %>%
    mutate( # summaries of the no PCRs and no + PCRS per blood sample
      N_PCRs_sample = length(CT),
      PCR_replicate = 1:length(CT),
      N_pos_sample = sum(CT<40)
    ) %>% ungroup() %>%
    arrange(ID_numeric, Day_frm_rand, MBREFID, MBGRPID) %>%
    group_by(ID_numeric) %>%
    mutate(Baseline_CT = mean(CT[which(Day_frm_rand <= 0)])) %>% 
    ungroup() %>%
    mutate(
      ARM = factor(ARM,levels = arms_select)
    ) %>%
    filter(PCR_replicate <= max_PCR_replicates) %>%
    group_by(ID_numeric, Day_frm_rand, MBREFID) %>%
    mutate( # recompute summaries after having filtered out measurements with more than 3 replicates
      N_PCRs_sample = length(CT),
      PCR_replicate = 1:length(CT),
      N_pos_sample = sum(CT<40)) %>%
    group_by(ID_numeric, Day_frm_rand) %>%
    mutate(
      N_PCRs = length(CT),
      N_pos = sum(CT<40))
  
  return(pcr_chagas_ss_input)
}



plot_individual_fits = function(data_list_stan,xx_CT_t,xx_CT, f_out='model_fits.pdf'){
  pdf(file = f_out, width = 12, height = 12)
  
  par(mfrow=c(3,3),las=1, family='serif',mar=c(2,2,2,1))
  for(id in 1:data_list_stan$N){
    ind = which(data_list_stan$id_ind==id)
    ts = data_list_stan$t_actual[ind]
    
    xlims = range(data_list_stan$t_actual)
    plot(NA, NA, xlim = xlims, 
         ylim = c(28, 50),
         xlab='',ylab='',panel.first=grid(),
         main=arms_all[data_list_stan$trt[id]])
    polygon(x = c(xlims+c(-100,100), rev(xlims)+c(100,-100)),
            y = c(40,40,100,100),border = NA, col = adjustcolor('grey',.2))
    abline(h=40)
    for(t in ind){
      for(j in 1:3){
        kmax=data_list_stan$k_replicates[t,j]
        if(kmax>0){
          for(k in 1:kmax){
            points(data_list_stan$t_actual[t], data_list_stan$CT_obs[t, j, k],
                   pch= 1 + as.numeric(data_list_stan$CT_obs[ind, j, k]==40)*16)
          }
        }
      }
    }
    lines(ts, xx_CT_t[ind], type='b',pch=16,lwd=1)
    lines(ts, xx_CT[ind], col='red',lwd=3)
  }
  dev.off()
}
