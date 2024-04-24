data {
  int<lower=0> N; // number of patients
  int<lower=1> K_arms; // number of randomised arms
  array[N] int<lower=1,upper=K_arms> trt_group; // trt arm for each patient
  int<lower=1> Kmax;  // max number of follow-up visits
  array[N] int<lower=0,upper=Kmax> y_pos; // number of visits with +PCR per patient
  array[N] int<lower=0,upper=Kmax> K_FUP; // number of visits per patient
}

parameters {
  vector<lower=0,upper=1>[K_arms] p_cured;
  vector<lower=0,upper=1>[2] PCR_sensitivity;
  real<lower=0,upper=1> PCR_specificity;
}


model {
  p_cured ~ beta(.5,.5);
  PCR_specificity ~ beta(500,1);
  PCR_sensitivity ~ beta(2,2);
  
  for(i in 1:N){
    real log_lik_given_cured;
    real log_lik_given_not_cured;
    real p_sens;
    if(trt_group[i]==1){
      p_sens = PCR_sensitivity[1];
    } else {
      p_sens = PCR_sensitivity[2];
    }
    log_lik_given_cured = binomial_lpmf(y_pos[i] | K_FUP[i], 0);
    log_lik_given_not_cured = binomial_lpmf(y_pos[i] | K_FUP[i], p_sens);
    target += 
    log_sum_exp(log_lik_given_cured+log(p_cured[trt_group[i]]), 
    log_lik_given_not_cured+ log1m(p_cured[trt_group[i]]));
  }
}

