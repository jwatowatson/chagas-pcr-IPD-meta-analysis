library(tidyverse)
library(cmdstanr)
# library(rstan)
library(posterior)
library(bayesplot)
library(gridExtra)
library(RColorBrewer)


file <- file.path("Stan_models/chagas_efficacy_model_CT_values_v4.stan")
mod <- cmdstan_model(file)

source('functions.R')
my_breaks = c(10, 30, 34, 36, 38, 39.99999,40)
my_breaks_legend = c('<30','30-34','34-36','36-38','38-40','>=40')
my_cols = c(brewer.pal(11, 'RdYlBu')[c(1,2,4,9,11)],brewer.pal(9, 'Pastel1')[9])


load('RData/pcr_chagas.RData')
pcr_chagas = pcr_chagas %>% arrange(ARM, USUBJID, Day_frm_rand, MBREFID, MBGRPID)%>%
  filter(!is.na(STUDYID), STUDYID=='BENDITA') %>%
  group_by(ARM) %>%
  mutate(
    BNZ_total_dose_allocated = case_when(
      ARM=='BNZ 300MG 8W' ~ 8*7*300,
      ARM=='BNZ 300MG 4W' ~ 4*7*300,
      ARM=='BNZ 300MG 2W' ~ 2*7*300,
      ARM=='BNZ 300MG/WK + E1224 300MG 8W' ~ 8*300,
      ARM=='BNZ 150MG 4W + E1224 300MG 8W' ~ 4*7*150,
      ARM=='BNZ 150MG 4W' ~ 4*7*150,
      T ~ 0
    ),
    BNZ_total_duration_allocated = case_when(
      ARM=='BNZ 300MG 8W' ~ 8*7,
      ARM=='BNZ 300MG 4W' ~ 4*7,
      ARM=='BNZ 300MG 2W' ~ 2*7,
      ARM=='BNZ 300MG/WK + E1224 300MG 8W' ~ 8*7,
      ARM=='BNZ 150MG 4W + E1224 300MG 8W' ~ 4*7,
      ARM=='BNZ 150MG 4W' ~ 4*7,
      T ~ 0
    )
  )

pcr_chagas_summary_t = pcr_chagas%>%
  distinct(STUDYID, USUBJID, Day_frm_rand, .keep_all = T)

pcr_chagas_summary = pcr_chagas%>%
  group_by(USUBJID) %>%
  mutate(FUP_time = max(Day_frm_rand) - min(Day_frm_rand)) %>%
  distinct(STUDYID, USUBJID, .keep_all = T) %>%
  group_by(USUBJID) %>%
  mutate(
    PP_pop = BNZ_total_dose >= 0.8*BNZ_total_dose_allocated,
    PP_pop2 = BNZ_total_days >= 0.8*BNZ_total_duration_allocated)

table(pcr_chagas_summary$ARM, pcr_chagas_summary$PP_pop)
table(pcr_chagas_summary$ARM, pcr_chagas_summary$PP_pop2)

PP_ids = pcr_chagas_summary$USUBJID[pcr_chagas_summary$PP_pop2]


table(pcr_chagas_summary$ARM, pcr_chagas_summary$BNZ_total_dose)
sum(pcr_chagas_summary$FUP_time)/365
table(pcr_chagas_summary$ARM)

table(pcr_chagas_summary$STUDYID, useNA = 'ifany')
unique(pcr_chagas_summary$USUBJID)

unique(pcr_chagas$ARM)

# pcr_chagas = pcr_chagas %>% filter(ARM=='PLACEBO')
arms_select = c('PLACEBO',
                'BNZ 300MG 8W',
                'BNZ 300MG/WK + E1224 300MG 8W',
                'BNZ 300MG 4W',
                "BNZ 150MG 4W",
                "BNZ 150MG 4W + E1224 300MG 8W",
                'BNZ 300MG 2W')#,
                # "E1224 400MG 4W",
                # "E1224 400MG 8W",
                # "E1224 200MG 8W")

pcr_chagas_ss_input = get_analysis_dataset(dat = pcr_chagas,
                                           arms_select = arms_select, 
                                           analysis_IDs = PP_ids,max_PCR_replicates = 3)

table(pcr_chagas_ss_input$STUDYID)
table(pcr_chagas_ss_input$PCR_replicate)



# Analysis sequence.
# Primary analysis is per protocol, whereby per protocol is defined as taking at least 80% of assigned doses
# Sensitivity analysis is ITT
# Separate analysis looking at placebo patients only
# does Student-t distribution help/better fit?


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


Nchains=1
fit_1 <- mod$sample(
  data = data_list_stan,
  seed = 123,
  chains = Nchains,
  parallel_chains = Nchains,
  iter_warmup=1000,
  iter_sampling = 2000, 
  init = 1,
  thin = 15,
  refresh = 500,
  output_dir='stan_out/',
  save_warmup=F
)
mcmc_trace(fit_1$draws(c(#'lambda_trt',
                        # 'p_1',
                         'sigma_PCR',
                         'sigma_sample',
                         'sigma_time',
                         'alpha_PCR',
                         'CT_sigma_pop_baseline',
                         'CT_pop_baseline')))
mcmc_intervals(fit_1$draws('eta_t'))
mcmc_intervals(fit_1$draws('eta_t_j'))
mcmc_intervals(fit_1$draws('CT_i'))

save(fit_1, file = 'stan_out/full_model.rds')

label_change = c("p_1[1]" ="BNZ 300mg 8wks",
                 "p_1[2]" ="BNZ 300 wkly + E1224",
                 "p_1[3]" ="BNZ 300mg 4wks",
                 "p_1[4]" ="BNZ 150mg 4wks",
                 "p_1[5]" ="BNZ 150mg 4wks + E1224",
                 "p_1[6]" ="BNZ 300mg 2wks",
                 "p_1[7]" ="E1224 HD 4W",
                 "p_1[8]" ="E1224 HD 8W",
                 "p_1[9]" ="E1224 LD 8W")

summary_p1 = summary(fit_1$draws('p_1'))
summary_p1$variable=label_change
summary_p1[, c('variable', 'mean','q5','q95')]

placebos_ct_i = paste('CT_i[',which(data_list_stan$trt==0),']', sep='')
xx=mcmc_intervals_data(fit_1$draws(placebos_ct_i))
par(mar=c(4,1,1,1),las=1)
plot(xx$m, 1:nrow(xx), xlim = c(min(xx$ll), max(xx$hh)),
     xlab='Ct value', panel.first = grid())
for(i in 1:nrow(xx)){
  lines(c(xx$ll[i], xx$hh[i]),c(i,i))
  lines(c(xx$l[i], xx$h[i]),c(i,i),lwd=3)
}



pdf('model_fits.pdf',width = 12, height = 12)


eta_t = colMeans(fit_1$draws('eta_t',format = 'matrix'))
ct_i = colMeans(fit_1$draws('CT_i',format = 'matrix'))

ps = (fit_1$draws('p_1',format = 'matrix'))
dim(ps)

par(las=1)
plot(data_list_stan$CT_blood_sample_mean, eta_t,
     ylab='Mean observed Ct value',
     xlab='Model predicted Ct value', 
     panel.first = grid())
lines(0:40, 0:40,lwd=3, col='red')


label_change2 = label_change
names(label_change2) =  gsub(pattern = 'p_1',
                             replacement = 'lambda_trt',x = names(label_change),fixed = T)

#mcmc_intervals(fit_1$draws(c('CT_i')))
mcmc_intervals(fit_1$draws('p_1'))+
  ggplot2::scale_y_discrete(labels =label_change)+
  ggtitle('Probability of cure')+xlim(0,1)+theme_minimal()

mcmc_intervals(fit_1$draws('lambda_trt'))+
  ggplot2::scale_y_discrete(labels =label_change2)+
  ggtitle('Change in steady state density if not cured (Ct scale)')+
  xlim(0,10)+theme_minimal()



par(mfrow=c(3,3),las=1, family='serif',mar=c(2,2,2,1))
for(id in 1:data_list_stan$N){
  ind = which(data_list_stan$id_ind==id)
  ts = data_list_stan$t_actual[ind]
  
  my_cols = c('purple','darkgreen','darkblue')
  xlims = range(data_list_stan$t_actual)
  plot(NA, NA, xlim = xlims, 
       ylim = c(28, 50),
       xlab='',ylab='',panel.first=grid(),
       main=arms_select[data_list_stan$trt[id]+1])
  polygon(x = c(xlims+c(-100,100), rev(xlims)+c(100,-100)),
          y = c(40,40,100,100),border = NA, col = adjustcolor('grey',.2))
  abline(h=40)
  #mtext(text = '',side = 1,line = 2)
  for(j in 1:3){
    for(k in 1:3){
      ct_vals=data_list_stan$CT_obs[ind, j, k]
      points(ts, ct_vals,col=my_cols[j],
             pch=j + as.numeric(ct_vals==40)*15)
    }
    
  }
  lines(ts, eta_t[ind]+ct_i[data_list_stan$id_ind[ind]], type='b',pch=16,lwd=1)
  lines(ts, rep(ct_i[id],length(ts)),lty=3,lwd=2)
  abline(v=ts[which(diff(data_list_stan$treatment[ind])>0)]+0.5)
}
dev.off()

