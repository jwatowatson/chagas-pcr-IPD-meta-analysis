library(tidyverse)
library(cmdstanr)
# library(rstan)
library(posterior)
library(bayesplot)
library(gridExtra)
library(RColorBrewer)


file <- file.path("Stan_models/chagas_efficacy_model_CT_values_v3.stan")
mod <- cmdstan_model(file)

source('functions.R')
my_breaks = c(10, 30, 34, 36, 38, 39.99999,40)
my_breaks_legend = c('<30','30-34','34-36','36-38','38-40','>=40')
my_cols = c(brewer.pal(11, 'RdYlBu')[c(1,2,4,9,11)],brewer.pal(9, 'Pastel1')[9])


load('RData/pcr_chagas.RData')
pcr_chagas = pcr_chagas %>% arrange(ARM, USUBJID, Day_frm_rand, MBREFID, MBGRPID)
#%>%
  #filter(STUDYID=='BENDITA')
unique(pcr_chagas$ARM)

arms_select = c('PLACEBO',
                'BNZ 300MG 8W',
                'BNZ 300MG/WK + E1224 300MG 8W',
                'BNZ 300MG 4W',
                "BNZ 150MG 4W",
                "BNZ 150MG 4W + E1224 300MG 8W",
                'BNZ 300MG 2W',
                "E1224 400MG 4W",
                "E1224 400MG 8W",
                "E1224 200MG 8W")

pcr_chagas_ss_input = pcr_chagas %>% 
  filter(ARM %in% arms_select,
         Day_frm_rand <=0 | #baseline
           Day_frm_rand > (EOT+14) | #follow-up
           ARM=='PLACEBO' # not treated
         ) %>%
  mutate(EOT = ifelse(ARM=='PLACEBO', 1000, EOT)) %>%
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
  mutate(Baseline_CT = mean(CT[which(VISIT=='Screening')])) %>% 
  ungroup() %>%
  mutate(
    ARM = factor(ARM,levels = arms_select)
    ) %>%
  filter(PCR_replicate <= 3)%>%
  group_by(ID_numeric, Day_frm_rand, MBREFID) %>%
  mutate( # recompute summaries after having filtered out measurements with more than 3 replicates
    N_PCRs_sample = length(CT),
    PCR_replicate = 1:length(CT),
    N_pos_sample = sum(CT<40)) %>%
  group_by(ID_numeric, Day_frm_rand) %>%
  mutate(
    N_PCRs = length(CT),
    N_pos = sum(CT<40))

table(pcr_chagas_ss_input$STUDYID)
table(pcr_chagas_ss_input$PCR_replicate)

data_list_stan = 
  make_stan_dataset_model_quant(pcr_chagas_ss_input = pcr_chagas_ss_input, 
                                p_1_beta_prior = .5)

table(n_pos = data_list_stan$y_pos, data_list_stan$K_total)
data_list_stan$Trt_max
data_list_stan$N
data_list_stan$Tmax
data_list_stan$treatment
data_list_stan$id_ind
data_list_stan$trt


Nchains=4
fit_1 <- mod$sample(
  data = data_list_stan,
  seed = 123,
  chains = Nchains,
  parallel_chains = Nchains,
  iter_warmup=2000,
  iter_sampling = 3000, 
  init = 1,
  thin = 15,
  refresh = 1000,
  output_dir='stan_out/',
  save_warmup=F
)

save(fit_1, file = 'stan_out/full_model.rds')
mcmc_trace(fit_1$draws(c('lambda_trt',
                         'p_1',
                         'sigma_PCR',
                         'sigma_sample',
                         'sigma_time',
                         'alpha_PCR',
                         'CT_sigma_pop_baseline','CT_pop_baseline')))

mcmc_intervals(fit_1$draws(c('CT_i')))
mcmc_intervals(fit_1$draws('p_1'))

ct_t = colMeans(fit_1$draws('CT_t',format = 'matrix'))
ct_i = colMeans(fit_1$draws('CT_i',format = 'matrix'))

plot(data_list_stan$CT_blood_sample_mean, ct_t,
     xlab='Mean observed CT value', ylab='Model predicted CT', 
     panel.first = grid())
lines(0:40, 0:40,lwd=3, col='red')

pdf('model_fits.pdf',width = 10, height = 10)
par(mfrow=c(2,2),las=1, family='serif',mar=c(2,2,2,1))
for(id in 1:data_list_stan$N){
  ind = which(data_list_stan$id_ind==id)
  ts = data_list_stan$t_actual[ind]
  
  my_cols = c('purple','pink','darkblue')
  plot(NA, NA, xlim = range(data_list_stan$t_actual), 
       ylim = c(28, 50),
       xlab='',ylab='',panel.first=grid(),
       main=arms_select[data_list_stan$trt[id]+1])
  #mtext(text = '',side = 1,line = 2)
  for(j in 1:3){
    for(k in 1:3){
      ct_vals=data_list_stan$CT_obs[ind, j, k]
      points(ts, ct_vals,col=my_cols[j],
             pch=j + as.numeric(ct_vals==40)*15)
    }
    lines(ts, ct_t[ind], type='b',pch=16,lwd=2)
  }
  lines(ts, rep(ct_i[id],length(ts)),lty=3,lwd=2)
  abline(v=ts[which(diff(data_list_stan$treatment[ind])>0)]+0.5)
}
dev.off()

