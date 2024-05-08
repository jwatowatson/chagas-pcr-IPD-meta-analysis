// modelling the CT values directly
// simplest model with binary treatment allocation: trt versus not trt
functions {
  real f(real mu, real alpha) {
    real sigma_het;
    real het_threshold=30.0;
    if(mu<het_threshold){
      sigma_het=1;
    } else {
      sigma_het = 1 + alpha*(mu-het_threshold);
    }
    return sigma_het;
  }
}


data {
  int<lower=0> N; // number of patients
  int<lower=N> Tmax; // max number of timepoints where blood was taken for all individuals
  int<lower=0,upper=3> Smax; // max number of blood samples
  int<lower=1,upper=3> Kmax; // max number of replicates per blood sample
  
  int<lower=0> Trt_max; // number of treatment arms
  
  array[N] int<lower=0,upper=Trt_max> trt; //0: placebo; 1..Trt_max indexes treatment given by individual
  //array[N] int<lower=0,upper=Tmax> n_follow_up; //number of timepoints after treatment
  array[Tmax] int<lower=0,upper=Smax> j_samples; //number of samples per timepoint
  array[Tmax,Smax] int<lower=0,upper=Smax> k_replicates; //number of samples per timepoint
  
  array[Tmax] int<lower=0,upper=1> treatment; // indexes treatment by timepoint (0 either placebo or before treatment)
  
  // indexing on the 1:Tmax dimension
  array[Tmax] int<lower=1,upper=N> id_ind;
  array[N] int<lower=1,upper=Tmax> ind_start;
  array[N] int<lower=1,upper=Tmax> ind_end;
  
  // CT value data
  array[Tmax,Smax,Kmax] real<lower=0,upper=40> CT_obs;
  
  // positive/negative data are used if cured (all should be negative)
  array[Tmax] int<lower=1,upper=9> K_total; // total number of PCRs done at timepoint t
  array[Tmax] int<lower=0,upper=9> y_pos; // number positive PCR at timepoint t
  
  // priors
  real<lower=0> p_1_beta_prior;
  //real<lower=0> nu_student_t;
}

transformed data {
  real CTmax = 50;
  real CTmin = 25;
}


parameters {
  real<lower=CTmin,upper=CTmax> CT_pop_baseline;         // mean population baseline value
  real<lower=0> CT_sigma_pop_baseline;                   // SD of population CT at baseline
  vector<lower=CTmin,upper=CTmax>[N] CT_i;               // individual baseline values
  vector<lower=CTmin,upper=CTmax>[Tmax] CT_t;            // mean CT value at timepoint t
  array[Tmax,Smax] real<lower=CTmin,upper=CTmax> CT_t_j; // mean CT value in sample j at timepoint t
  
  // treatment effects
  real<lower=-10,upper=10> lambda_trt;        // treatment effect parameterised as a shift in mean CT value
  array[Trt_max] real<lower=0,upper=1> p_1; // probability of full cure
  
  // measurement error
  real<lower=0> sigma_PCR;
  real<lower=0> sigma_sample;
  real<lower=0> sigma_time;
  real<lower=0,upper=2> alpha_PCR;            // heteroskedasticity parameter for PCR measurement
}


model {
  
  //prior
  p_1 ~ beta(p_1_beta_prior,p_1_beta_prior);
  lambda_trt ~ normal(0.0, 2.0) T[-10,10];
  sigma_PCR ~ normal(0.5, 0.25) T[0,];
  alpha_PCR ~ normal(0.5, 0.25) T[0,];
  sigma_sample ~ normal(2, 0.5) T[0,];
  sigma_time ~ normal(2, 0.5) T[0,];
  
  // hierachical model of the baseline values
  CT_pop_baseline ~ normal(40, 2) T[CTmin, CTmax];
  CT_sigma_pop_baseline ~ normal(4, 1) T[0,];
  CT_i ~ normal(CT_pop_baseline, CT_sigma_pop_baseline) T[CTmin,CTmax];
  
  for(t in 1:Tmax){
    if(treatment[t]==0){
      CT_t[t] ~ normal(CT_i[id_ind[t]], sigma_time) T[CTmin,CTmax];
    } else {
      CT_t[t] ~ normal(CT_i[id_ind[t]]+lambda_trt, sigma_time) T[CTmin,CTmax];
    }
    for(j in 1:j_samples[t]){
      CT_t_j[t,j] ~ normal(CT_t[t], sigma_sample) T[CTmin,CTmax];
    }
  }
  
  //likelihood
  for(i in 1:N){
    // not treated case: assume constant mean value over time
    if(trt[i]==0){ 
      for(t in ind_start[i]:ind_end[i]){
        for(j in 1:j_samples[t]){
          real sigma_t_j = sigma_PCR*f(CT_t_j[t,j], alpha_PCR);
          for(k in 1:k_replicates[t,j]){
            if(CT_obs[t,j,k]<40){
              //observed
              target += student_t_lpdf(CT_obs[t,j,k] | 10, CT_t_j[t,j], sigma_t_j);
            } else {
              //censored at 40
              target += student_t_lccdf(40 | 10, CT_t_j[t,j], sigma_t_j);
            }
          }
        }
      }
    } else { // trt case: (i) cured or (ii) not cured: function of expected CT value
      
      real log_lik_given_cured=0;
      real log_lik_given_not_cured=0;
      // iterate over timepoints
      for(t in ind_start[i]:ind_end[i]){
        //**** if cured *****
        // before treatment
        if(treatment[t]==0){
          for(j in 1:j_samples[t]){
            real sigma_t_j = sigma_PCR*f(CT_t_j[t,j], alpha_PCR);
            for(k in 1:k_replicates[t,j]){
              if(CT_obs[t,j,k]<40){
                //observed
                log_lik_given_cured += student_t_lpdf(CT_obs[t,j,k] | 10, CT_t_j[t,j], sigma_t_j);
              } else {
                //censored at 40
                log_lik_given_cured += student_t_lccdf(40 | 10, CT_t_j[t,j], sigma_t_j);
              }
            }
          }
        } else {
          // after treatment
          log_lik_given_cured += binomial_lpmf(y_pos[t] | K_total[t], 0);
        }

        //**** if NOT cured *****
        for(j in 1:j_samples[t]){
          real sigma_t_j = sigma_PCR*f(CT_t_j[t,j], alpha_PCR);
          for(k in 1:k_replicates[t,j]){
            if(CT_obs[t,j,k]<40){
              //observed
              log_lik_given_not_cured += student_t_lpdf(CT_obs[t,j,k] | 10, CT_t_j[t,j], sigma_t_j);
            } else {
              //censored at 40
              log_lik_given_not_cured += student_t_lccdf(40 | 10, CT_t_j[t,j], sigma_t_j);
            }
          }
        }
      }
      // marginalise over the cured/not cured outcomes
      target +=
      log_sum_exp(log_lik_given_cured+log(p_1[trt[i]]),
      log_lik_given_not_cured+log1m(p_1[trt[i]]));
    }
  }
}
