data {
  int<lower=0> n_id; // number of patients
  int<lower=n_id> n_samples; // number of blood samples taken
  int<lower=1> K_arms; // number of randomised arms
  array[n_id] int<lower=1,upper=K_arms> trt_group; // trt arm for each patient
  
  int K_cov;
  matrix[n_id, K_cov] X;
  
  int<lower=1> Kmax; // max number of replicates per blood sample
  array[n_id] int<lower=1,upper=n_samples> ind_start;
  array[n_id] int<lower=1,upper=n_samples> ind_end;
  array[n_samples] int<lower=0,upper=Kmax> y_pos; // number positive PCR per blood sample
  array[n_samples] int<lower=1,upper=Kmax> K_FUP; // number of PCRs replicates per blood sample
  
  // tuning variable - specificity of the PCR at given cutoff threshold
  real<lower=0,upper=1> PCR_specificity;
}


parameters {
  vector[K_arms] logit_p_1; // probability of full cure
  vector[K_arms] logit_p_2; // probability of sampling parasites if not fully cured
  real<lower=0,upper=1> q; // sensitivity of individual PCR
  vector[K_cov] beta; //covariate effect on full cure (p_1)
}

transformed parameters {
  vector[n_id] p_1_id;
  vector[K_arms] p_1;
  vector[K_arms] p_2;
  vector[n_id] x_pred = X*beta;
  for(i in 1:n_id){
    p_1_id[i] = inv_logit(x_pred[i] + logit_p_1[trt_group[i]]);
  }
  p_1 = inv_logit(logit_p_1);
  p_2 = inv_logit(logit_p_2);
}


model {
  logit_p_1 ~ normal(0,1);
  logit_p_2 ~ normal(0,1);
  q ~ beta(2,2);
  beta ~ normal(0,1);
  
  for(i in 1:n_id){
    real log_lik_given_cured=0;
    real log_lik_given_not_cured=0;
    
    // take all PCRs as independent (might be some batch dependencies...)
    log_lik_given_cured = binomial_lpmf(sum(y_pos[ind_start[i]:ind_end[i]]) | sum(K_FUP[ind_start[i]:ind_end[i]]), 1-PCR_specificity);

    for(j in ind_start[i]:ind_end[i]){
      log_lik_given_not_cured += 
      log_sum_exp(log(p_2[trt_group[i]]) + binomial_lpmf(y_pos[j] | K_FUP[j], q), 
      log1m(p_2[trt_group[i]]) + binomial_lpmf(y_pos[j] | K_FUP[j], 1-PCR_specificity));
    }
    target += 
      log_sum_exp(log_lik_given_cured+log(p_1_id[i]), 
      log_lik_given_not_cured+log1m(p_1_id[i]));
  }
}

