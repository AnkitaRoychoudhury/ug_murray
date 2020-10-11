data {
  int<lower=0> N;
  real k[N];
  int<lower=0> N_ppc;
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

model {

    //Priors
  alpha ~ lognormal(0, 2);
  beta ~ lognormal(0,3);
  //Likelihood
  
  //k ~ normal(0, 1);
  
  k ~ gamma(alpha, beta);
}

generated quantities{

    real k_ppc[N_ppc];
    real log_lik[N];


    for (i in 1:N_ppc){
        k_ppc[i] = gamma_rng(alpha, beta);
        }
        
    for (i in 1:N) {
        log_lik[i] = gamma_lpdf(k[i] | alpha, beta);
        }
}
