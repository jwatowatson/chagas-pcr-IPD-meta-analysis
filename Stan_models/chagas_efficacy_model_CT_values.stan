// modelling the CT values directly
// simplest model with binary treatment allocation: treated versus not treated
functions {
  real f(real mu, real alpha) {
    real sigma_het;
    if(mu<30){
      sigma_het=0;
    } else {
      sigma_het = alpha*(mu-30);
    }
    return sigma_het;
  }
}


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
  real<lower=0, upper=45> CT_pop_baseline;
  real<lower=0> CT_sigma_pop_baseline;
  vector<lower=0, upper=45>[n_id] CT_i_baseline;
  vector<lower=0, upper=50>[n_samples] CT_blood_sample; // theoretical mean CT value in sample j
  real lambda_trt; // treatment effect parameterised as a shift in mean CT value
  real<lower=0,upper=1> p_1; // probability of full cure
  real<lower=0> sigma_PCR;
  real<lower=0> sigma_sample;
  real<lower=0> alpha_PCR;
  real<lower=0> alpha_sample;
}

transformed parameters {
  vector[n_samples] CT_timepoint; // mean CT value at each timepoint in those not cured
  for(j in 1:n_samples){
    CT_timepoint[j] = CT_i_baseline[id_ind[j]] + lambda_trt * treatment[j];
  }
}

model {
  
  //prior
  p_1 ~ beta(p_1_beta_prior,p_1_beta_prior);
  lambda_trt ~ normal(0.0 , 2);
  sigma_PCR ~ normal(0.5, 0.5) T[0,];
  alpha_sample ~ normal(0, 1.0) T[0,];
  alpha_PCR ~ normal(0, 1.0) T[0,];

  sigma_sample ~ normal(0, 2.0) T[0,];
  
  CT_pop_baseline ~ normal(35.0, 3.0);
  CT_sigma_pop_baseline ~ normal(2.0, 2.0) T[0,];
  CT_i_baseline ~ normal(CT_pop_baseline, CT_sigma_pop_baseline);
  
  for(i in 1:n_samples){
    CT_blood_sample[i] ~ normal(CT_timepoint[i], sigma_sample+f(CT_timepoint[i],alpha_sample));
  }
  
  //likelihood
  for(i in 1:n_id){
    
    if(treated[i]==0){ // not treated case: simply a function of expected CT value
      for(j in ind_start[i]:ind_end[i]){
        for(k in 1:Kmax){
          if(CT_obs[j,k]<40){//observed
          target += normal_lpdf(CT_obs[j,k] | CT_blood_sample[j], sigma_PCR + sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
          } else {//censored at 40
          target += normal_lccdf(40 | CT_blood_sample[j], sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
          }
        }
      }
      
    } else 
    {// treated case: (i) cured or (ii) not cured: function of expected CT value
      real log_lik_given_cured=0;
      real log_lik_given_not_cured=0;
      for(j in ind_start[i]:ind_end[i]){
        // if cured
        if(treatment[j]==0){
          // before treatment
          for(k in 1:Kmax){
            if(CT_obs[j,k]<40){//observed
            log_lik_given_cured += normal_lpdf(CT_obs[j,k] | CT_blood_sample[j], sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
            } else {//censored at 40
            log_lik_given_cured += normal_lccdf(40 | CT_blood_sample[j], sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
            }
          }
        } else {
          // after treatment
          log_lik_given_cured += binomial_lpmf(y_pos[j] | K_FUP[j], 0);
        }
        // if not cured
        for(k in 1:Kmax){
          if(CT_obs[j,k]<40){//observed
          log_lik_given_not_cured += normal_lpdf(CT_obs[j,k] | CT_blood_sample[j], sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
          } else {//censored at 40
          log_lik_given_not_cured += normal_lccdf(40 | CT_blood_sample[j], sigma_PCR+f(CT_blood_sample[j],alpha_PCR));
          }
        }
      }
      // marginalise over the cured/not cured outcomes
      target += 
      log_sum_exp(log_lik_given_cured+log(p_1), 
      log_lik_given_not_cured+log1m(p_1));
    }
  }
}

