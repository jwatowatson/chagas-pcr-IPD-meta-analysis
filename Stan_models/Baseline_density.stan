// this estimates baseline densities under a Poisson model
data {
  int<lower=1> N; // number of individuals
  int<lower=1> I; // number of blood samples
  int<lower=1> J; // number of repeats
  
  array[N,I,J] real<lower=0,upper=40> Y; // observed blood samples - this 40-CT (inverted)
  array[N,I,J] int<lower=0,upper=1> Y_40; // censored yes=1/no=0
  
  real<lower=0> vol_ml_sample;
  int max_parasites;
}


parameters {
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=40> beta;
  vector[N] log_lambda; // log density
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mean_parasite_dens;
  mean_parasite_dens = vol_ml_sample*exp(log_lambda);
}

model {
  
  beta ~ normal(38, 1);
  alpha ~ beta(95,5);
  log_lambda ~ normal(0,2); // mean is 1 parasite per ml
  sigma ~ normal(1,1);
  
  for(n in 1:N){
    
    for(i in 1:I){
      
      vector[max_parasites+1] log_lik_n_parasites;
      vector[J] log_lik_repeat;
      
      // special case of zero parasites - assume 100% specificity of PCR
      log_lik_n_parasites[1] = poisson_log_lpmf(0 | log(vol_ml_sample) + log_lambda[n]);
      for(j in 1:J){
        if(Y_40[n,i,j]==1){
          log_lik_repeat[j] = 0;
        } else {
          log_lik_repeat[j] = log(0); // assumes 100% specificity if no parasites
        }
      }
      log_lik_n_parasites[1] += sum(log_lik_repeat);
      
      // iterate over latent number of parasites in the ith blood draw
      for(k in 2:(max_parasites+1)){
        log_lik_n_parasites[k] = poisson_log_lpmf(k-1 | log(vol_ml_sample) + log_lambda[n]);
        for(j in 1:J){
          if(Y_40[n,i,j]==1){
            log_lik_repeat[j]=normal_lccdf(40 | -log2(k-1)+beta, sigma);
          } else {
            log_lik_repeat[j]=normal_lpdf(Y[n,i,j] | -log2(k-1)+beta, sigma);
          }
        }
        log_lik_n_parasites[k] += sum(log_lik_repeat);
      }
      
      target += log_sum_exp(log_lik_n_parasites);
      
    }
  }
}

generated quantities {
  array[N] int n_parasites;
  for(n in 1:N){
    n_parasites = poisson_rng(vol_ml_sample*exp(log_lambda));
  }
}

