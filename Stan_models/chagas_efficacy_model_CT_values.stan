// modelling the CT values directly
// simplest model with just binary treatment allocation

data {
  int<lower=0> n_id; // number of patients
  int<lower=n_id> n_samples; // number of blood samples taken
  array[n_id] int<lower=0,upper=1> treated; //0: not treated; 1: treated
  array[n_samples] int<lower=0,upper=1> treatment; //0: before treatment; 1: after treatment
  int<lower=1> Kmax; // max number of replicates per blood sample
  
  // indexing
  array[n_samples] int<lower=1,upper=n_id> id_ind;
  array[n_id] int<lower=1,upper=n_samples> ind_start;
  array[n_id] int<lower=1,upper=n_samples> ind_end;
  
  // CT value data are used if not cured
  matrix<lower=0,upper=40>[n_samples, Kmax] CT_obs;
  // positive/negative data are used if cured (all should be negative)
  array[n_samples] int<lower=0,upper=Kmax> y_pos; // number positive PCR per blood sample
  array[n_samples] int<lower=1,upper=Kmax> K_FUP; // number of PCRs replicates per blood sample
  
  // priors
  real<lower=0> p_1_beta_prior;
}


parameters {
  real CT_pop_baseline;
  real<lower=0> CT_sigma_pop_baseline;
  vector[n_id] CT_i_baseline;
  vector[n_samples] CT_sample; // theoretical mean CT value in sample j
  real lambda_trt; // treatment effect parameterised as a shift in mean CT value
  vector[n_samples] CT_i_sample;
  real<lower=0,upper=1> p_1; // probability of full cure
  real<lower=0> sigma_PCR;
  real<lower=0> sigma_sample;
}

transformed parameters {
  vector[n_samples] CT_sample_i_theoretical; // theoretical mean CT value at each timepoint
  for(j in 1:n_samples){
    CT_sample_i_theoretical[j] = CT_i_baseline[id_ind[j]] + lambda_trt * treatment[j];
  }
}

model {
  
  p_1 ~ beta(p_1_beta_prior,p_1_beta_prior);
  CT_pop_baseline ~ normal(35, 3);
  CT_sigma_pop_baseline ~ normal(2, 2);
  CT_i_baseline ~ normal(CT_pop_baseline, CT_sigma_pop_baseline);
  sigma_sample ~ normal(2,2);
  CT_sample ~ normal(CT_sample_i_theoretical, sigma_sample);
  
  for(i in 1:n_id){
    
    if(treated[i]==0){
      for(j in ind_start[i]:ind_end[i]){
        target += normal_lpdf(CT_obs[j,] | CT_sample[j], sigma_PCR);
      }
      
    } else {
      real log_lik_given_cured=0;
      real log_lik_given_not_cured=0;
      for(j in ind_start[i]:ind_end[i]){
        // if cured
        if(treatment[j]==0){
          // before treatment
          log_lik_given_cured += normal_lpdf(CT_obs[j,] | CT_sample[j], sigma_PCR);
        } else {
          // after treatment
          log_lik_given_cured += binomial_lpmf(y_pos[j] | K_FUP[j], 0);
        }
        
        // if not cured
        log_lik_given_not_cured += normal_lpdf(CT_obs[j,] | CT_sample[j], sigma_PCR);
      }
      target += 
      log_sum_exp(log_lik_given_cured+log(p_1*treatment[i]), 
      log_lik_given_not_cured+log1m(p_1*treatment[i]));
    }
  }
}

