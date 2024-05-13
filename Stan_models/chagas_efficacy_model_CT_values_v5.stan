// modelling the CT values directly
functions {
  real f(real mu, real alpha) {
    real sigma_het = -1;
    real het_threshold_lower = 30.0;
    real het_threshold_upper = 45;
    
    if(mu < het_threshold_lower){
      sigma_het = 1;
    } else if (mu > het_threshold_upper) {
      sigma_het = 1 + alpha*(het_threshold_upper-het_threshold_lower);
    } else {
      sigma_het = 1 + alpha*(mu-het_threshold_lower);
    }
    
    return sigma_het;
  }
  
  real compute_log_lik(real x, real mu, real sigma){
    real log_lik;
    if(x<40){  //observed
      log_lik = student_t_lpdf(x | 10, mu, sigma);
    } else {
      log_lik = student_t_lccdf(40 | 10, mu, sigma);
    }
    return(log_lik);
  }
}


data {
  int<lower=0> N; // number of patients
  int<lower=N> Tmax; // max number of timepoints where blood was taken for all individuals
  int<lower=0,upper=3> Smax; // max number of blood samples
  int<lower=1,upper=3> Kmax; // max number of replicates per blood sample
  
  int<lower=0> Trt_max; // number of treatment arms
  
  array[N] int<lower=1,upper=Trt_max> trt; // 1..Trt_max indexes allocated treatment by individual
  array[Tmax] int<lower=0,upper=Smax> j_samples; // number of samples per timepoint
  array[Tmax,Smax] int<lower=0,upper=Kmax> k_replicates; // number of samples per timepoint
  
  array[Tmax] int<lower=0,upper=1> EOT; // 0 before EOT, 1 after
  
  // indexing on the 1:Tmax dimension
  array[Tmax] int<lower=1,upper=N> id_ind;
  array[N] int<lower=1,upper=Tmax> ind_start;
  array[N] int<lower=1,upper=Tmax> ind_end;
  
  // CT value data
  array[Tmax,Smax,Kmax] real<lower=0,upper=40> CT_obs;
  
  // positive/negative data are used if cured (all should be negative)
  array[Tmax] int<lower=1,upper=9> K_total; // total number of PCRs done at timepoint t
  array[Tmax] int<lower=0,upper=9> y_pos;   // number positive PCR at timepoint t
  
  // priors
  real<lower=0> p_1_beta_prior;
  real lambda_max;
}

transformed data {
  real CTmax = 50;
  real CTmin = 25;
}


parameters {
  real<lower=CTmin,upper=CTmax> CT_pop_baseline;  // mean population baseline value
  real<lower=0> CT_sigma_pop_baseline;            // SD of population CT at baseline
  vector<lower=CTmin,upper=CTmax>[N] CT_i;        // individual baseline values
  
  // random effects terms
  array[Tmax] real eta_t;               // CT random effects by timepoint
  vector[Smax] eta_t_j[Tmax];           // CT random effect for sample j at timepoint t

  // treatment effects
  real<lower=0> sigma_trt_all;
  array[Trt_max] real<lower=-10,upper=lambda_max> lambda_trt;     // treatment effect parameterised as a shift in mean CT value
  array[Trt_max] real<lower=0,upper=1> p_1;               // probability of cure

  // measurement error
  real<lower=0> sigma_PCR;
  positive_ordered[2] sigma;     // sigma[1] is the sample variance; sigma[2] is the temporal variance
  real<lower=0> alpha_PCR;       // heteroskedasticity parameter for PCR measurement
}

transformed parameters {
  
  array[Tmax,Smax] real CT_t_j; // mean CT value in sample j at timepoint t
  array[Tmax] real CT_i_t;      // mean CT value at timepoint t without random effect
  array[Tmax] real CT_t;        // mean CT value at timepoint t with random effect
  vector[N] CT_i_EOT;        // individual values after EOT

  for(i in 1:N){
    CT_i_EOT[i] = CT_i[i]+lambda_trt[trt[i]];
  }
  for(t in 1:Tmax){
    
    if(EOT[t]==0) {
      CT_i_t[t] = CT_i[id_ind[t]];     
    } else {
      CT_i_t[t] = CT_i_EOT[id_ind[t]];
    }
    // add temporal random effect
    CT_t[t] = CT_i_t[t] + eta_t[t];
    // add sample random effect
    for(j in 1:Smax){
      CT_t_j[t,j] = CT_t[t] + eta_t_j[t][j];
    }
  }
}

model {
  
  //prior
  p_1 ~ beta(p_1_beta_prior,p_1_beta_prior);
  sigma_trt_all ~ normal(0,3) T[0,];
  lambda_trt ~ normal(0, sigma_trt_all);
  sigma_PCR ~ normal(0, 0.5) T[0,];
  alpha_PCR ~ normal(0, 0.5) T[0,];
  sigma[1] ~ normal(1,1) T[0,];
  sigma[2] ~ normal(2,1) T[0,];
  
  // hierachical model of the baseline values
  CT_pop_baseline ~ normal(40, 5) T[CTmin, CTmax];
  CT_sigma_pop_baseline ~ normal(0, 4) T[0,];
  
  // Individual steady states - contrained values
  CT_i ~ normal(CT_pop_baseline, CT_sigma_pop_baseline) T[CTmin,CTmax];
  
  // random effects
  eta_t ~ normal(0, sigma[2]);
  for(t in 1:Tmax){
    eta_t_j[t] ~ normal(0,sigma[1]); //multi_normal_cholesky(zeros2, diag_pre_multiply(sigma_u,Omega_chol));
  }
 
  //likelihood
  for(i in 1:N){
    real log_lik_given_cured=0;
    real log_lik_given_not_cured=0;
    // iterate over timepoints
    for(t in ind_start[i]:ind_end[i]){
      // **** if cured scenario *****
      if(EOT[t]==0){
        for(j in 1:j_samples[t]){
          real sigma_t_j = sigma_PCR*f(CT_t_j[t,j], alpha_PCR);
          for(k in 1:k_replicates[t,j]){
            log_lik_given_cured += compute_log_lik(CT_obs[t,j,k], CT_t_j[t,j], sigma_t_j);
          }
        }
      } else {
        // after EOT
        log_lik_given_cured += binomial_lpmf(y_pos[t] | K_total[t], 0);
      }
      
      //**** if NOT cured *****
      for(j in 1:j_samples[t]){
        real sigma_t_j = sigma_PCR*f(CT_t_j[t,j], alpha_PCR);
        for(k in 1:k_replicates[t,j]){
          log_lik_given_not_cured += compute_log_lik(CT_obs[t,j,k], CT_t_j[t,j], sigma_t_j);
        }
      }
    }
    // marginalise over the cured/not cured outcomes
    target +=
    log_sum_exp(log_lik_given_cured+log(p_1[trt[i]]),
    log_lik_given_not_cured+log1m(p_1[trt[i]]));
  }
}
