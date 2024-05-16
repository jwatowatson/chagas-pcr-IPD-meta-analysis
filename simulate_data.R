require(tidyverse)
require(foreach)

f_het = function(mu, min_val, max_val, alpha_pcr){
  sigma_het = 1
  if(mu>min_val & mu<max_val){
    sigma_het=1+alpha_pcr*(mu-min_val)
  }
  if(mu>max_val){
    sigma_het=1+alpha_pcr*(max_val-min_val)
  }
  return(sigma_het)
}


simulate_individual = 
  function(id=1, 
           N_visits_no_trt=2,
           N_visits_EOT=5, 
           p_cure=.8,
           N_blood_samples=2,
           N_PCR_replicates=3,
           xi=40, 
           lambda=0,
           sigma_t=2,
           sigma_s=1, 
           sigma_PCR=0.25, 
           alpha_PCR=1,
           ARM=1){
    
    if(length(id) != 1) stop('this function only simulates one individual')
    
    dat_SDTM = expand.grid(USUBJID=id, 
                           VISIT=1:(N_visits_no_trt+N_visits_EOT), 
                           MBREFID=1:N_blood_samples, 
                           MBGRPID=1:N_PCR_replicates,
                           xi=xi,
                           EOT = N_visits_no_trt,
                           STUDYID='SIM',
                           ARM=ARM)
    
    dat_SDTM = dat_SDTM %>% 
      mutate(Cured = rbinom(1,1,p_cure),
             Day_frm_rand=VISIT,
             VISIT_numeric=VISIT,
             VISIT_trans=VISIT) %>%
      arrange(USUBJID, VISIT, MBREFID, MBGRPID) %>%
      group_by(VISIT) %>%
      mutate(
        eta_t = rnorm(n = 1, 0, sigma_t),
        CT_t = xi + eta_t,
        CT_t = ifelse(VISIT>N_visits_no_trt, CT_t+lambda, CT_t)) %>%
      group_by(VISIT, MBREFID) %>%
      mutate(
        eta_t_j = rnorm(n=1, mean=0, sd=sigma_s),
        CT_t_j = CT_t + eta_t_j) %>% 
      group_by(VISIT, MBREFID, MBGRPID) %>%
      mutate(
        eta_PCR = rnorm(n=1, mean=0, sd = sigma_PCR*f_het(CT_t_j, 30, 45, alpha_PCR)),
        CT_t_j_k = CT_t_j + eta_PCR,
        CT = min(40, CT_t_j_k),
        CT = ifelse(Cured==1 & VISIT>N_visits_no_trt, 40, CT)
      ) %>%
      group_by(VISIT, MBREFID) %>%
      mutate(mean_CT = mean(CT),
             sd_CT = sd(CT))
    
    return(dat_SDTM)
    
  }

# View(simulate_individual())

simulate_pop = function(N=1,
                        N_visits_no_trt=2,
                        N_visits_EOT=5,
                        N_blood_samples=3,
                        N_PCR_replicates=3,
                        thetas,
                        p_index=1)
{
  ind_posterior = sample(1:length(thetas$CT_pop_baseline), size = N, replace = T)
  dat_sim = foreach(id = 1:N, .combine = rbind) %do% {
    j = ind_posterior[id]
    simulate_individual(id = id,
                        N_visits_no_trt = N_visits_no_trt, 
                        N_visits_EOT = N_visits_EOT,
                        p_cure = thetas$p_1[j,p_index ],
                        N_blood_samples = N_blood_samples,
                        N_PCR_replicates = N_PCR_replicates,
                        xi=rnorm(1, mean = thetas$CT_pop_baseline[j], 
                                 sd = thetas$CT_sigma_pop_baseline[j]), 
                        lambda = thetas$lambda_trt[j, p_index],
                        sigma_t = thetas$sigma[j,2],
                        sigma_s = thetas$sigma[j,1], 
                        sigma_PCR = thetas$sigma_PCR[j], 
                        alpha_PCR = thetas$alpha_PCR[j],
                        ARM=p_index)
  }
  return(dat_sim)
}

source('functions.R')
fname='stan_out/all_data_PP.rds'
load(fname)
thetas = rstan::extract(fit_1, pars=c('CT_pop_baseline','CT_sigma_pop_baseline',
                                      'sigma','sigma_PCR','alpha_PCR','p_1',
                                      'lambda_trt'))

sim_dat = simulate_pop(N = 30, thetas = thetas,p_index = 1)
pcr_chagas = sim_dat %>% 
  group_by(USUBJID, Day_frm_rand, MBREFID) %>%
  mutate(
    N_technical_replicates = 1:n()) %>%
  ungroup() %>%
  filter(N_technical_replicates<=3) %>%
  mutate(N_pos_sample = sum(CT<40)) %>%
  group_by(USUBJID,Day_frm_rand) %>% 
  arrange(USUBJID, Day_frm_rand, MBREFID, MBGRPID) %>%
  mutate(PCR_number = 1:n()) %>%
  ungroup() %>% 
  group_by(USUBJID,Day_frm_rand) %>%
  mutate(percent_pos = sum(CT<40)/length(CT),
         Any_Pos_40 = any(CT<40),
         N_PCRs = max(PCR_number),
         N_pos = sum(CT<40)) %>%
  group_by(USUBJID) %>%
  mutate(
    Percent_pos_baseline = mean(CT[Day_frm_rand<=1] < 40),
    Mean_CT_inv_baseline = 40-mean(CT[Day_frm_rand<=1]),
    N_FUP_VISITS = length(unique(Day_frm_rand[which(Day_frm_rand>EOT)])),
    Day_last_Visit = max(Day_frm_rand)) %>%
  ungroup() %>%
  filter(PCR_number<=9) %>%
  arrange(Mean_CT_inv_baseline, Percent_pos_baseline,
          USUBJID, VISIT_numeric, PCR_number)

my_breaks = c(10, 30, 34, 36, 38, 39.99999,40)
my_breaks_legend = c('<30','30-34','34-36','36-38','38-40','>=40')
my_cols = c(brewer.pal(11, 'RdYlBu')[c(1,2,4,9,11)],brewer.pal(9, 'Pastel1')[9])

pcr_mat = make_matrix_pcr(pcr_dat = pcr_chagas)
par(mar=c(3,8,1,8))
xx2 = pcr_chagas %>% ungroup %>% distinct(USUBJID,.keep_all = T)
ind_y = which(!duplicated(xx2$ARM))
plot_pcr_matrix(pcr_mat = pcr_mat,
                my_breaks = my_breaks, 
                my_break_legend = my_breaks_legend,
                my_cols = my_cols,
                arm_labels = unique(pcr_chagas$ARM),
                h_lines_ind = ind_y[-1]-0.5,
                y_lab_ind = ind_y + diff(c(ind_y, nrow(pcr_mat)+1))/2)
