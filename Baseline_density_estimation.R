library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(gridExtra)

file <- file.path("Baseline_density.stan")
mod <- cmdstan_model(file)

load('RData/pcr_chagas.RData')

dat = pcr_chagas %>% filter(ARM=='PLACEBO',VISIT=='Screening', USUBJID<100)

N = length(unique(dat$ID))
Y = Y_40 = array(dim = c(N, 3, 3))
dat$ID_ind = as.numeric(as.factor(dat$ID))
for(n in 1:N){
  for(i in 1:3){
    for(j in 1:3){
      ind = which(dat$ID_ind==n & dat$MBGRPID==i & dat$MBREFID==paste('Sample',j))
      if(length(ind) != 1){
        print(n)
      }
      Y[n,i,j]=40-dat$CT[ind]
      Y_40[n,i,j] = as.numeric(dat$CT[ind]==40)
    }
  }
}

data_list = list(
  N = N,
  I = 3,
  J = 3,
  Y = Y,
  Y_40 = Y_40,
  vol_ml_sample=5,
  max_parasites=1000
)


fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,iter_sampling = 2000,
  refresh = 500 # print update every 500 iters
)


mcmc_trace(fit$draws(c('alpha','beta','sigma')))
mcmc_intervals(fit$draws('log_lambda'))

par(las=1)
plot(log2(1:128), -(1/0.99)*log2(1:128)+39.5, type='l',
     xaxt='n',panel.first=grid(),ylim =c(30,40),
     xlab='parasites/ml', ylab='CT')
axis(1, at= 0:7, labels = 2^(0:7))

xx = apply(fit$draws('log_lambda'), 3, mean)
dat = dat %>% distinct(ID, .keep_all = T)
plot(dat$Mean_CT, log2(5*exp(xx)), panel.first=grid(),yaxt='n',
     xlab='Mean observed CT (9 values)',
     ylab='Estimated number of parasites in 5 mls')
axis(2, at= -1:6, labels = 2^(-1:6))

hist(xx, breaks = 10)
mu_dens = fit$draws('mean_parasite_dens')[,,12]
hist(mu_dens, breaks = 100)
xx = fit$draws('n_parasites')
xx = apply(xx, c(1,3), mean)
xx = rowMeans(fit$draws('mean_parasite_dens')[,,1, drop=T])
