data {
  int<lower=0> N;
  real k[N];
  int N_ppc;
}

parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
}

model {

    //Priors
  mu ~ lognormal(0, 10);
  sigma ~ lognormal(0,10);
  //Likelihood
  
  //k ~ normal(0, 1);
  
  k ~ normal(mu, sigma);
}

generated quantities{

    real k_ppc[N_ppc];
    real log_lik[N];
    
    for (i in 1:N_ppc){
        k_ppc[i] = normal_rng(mu, sigma);
        }
        
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(k[i] | mu, sigma);
        }
}

