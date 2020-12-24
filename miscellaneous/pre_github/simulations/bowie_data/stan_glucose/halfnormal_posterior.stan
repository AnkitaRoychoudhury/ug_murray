data {
  int<lower=0> N;
  real k[N];
  int<lower=0> N_ppc;
}

parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
}

model {
   
    //Priors
  mu ~ lognormal(0, 20);
  sigma ~ lognormal(0,30);
  //Likelihood
 //real<lower=mu> k;
  
  
  
  k ~ normal(mu, sigma);
}

generated quantities{

    real k_ppc[N_ppc];
    //real<lower=mu> k;
    for (i in 1:N_ppc){
        k_ppc[i] = mu + abs(normal_rng(0, sigma));


        }
}
