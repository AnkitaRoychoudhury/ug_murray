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
  mu ~ lognormal(0, 2);
  sigma ~ lognormal(2,3);
  //Likelihood
  
  //k ~ normal(0, 1);
  
  k ~ normal(mu, sigma);
}

generated quantities{

    real k_ppc[N_ppc];
    for (i in 1:N_ppc){
        k_ppc[i] = normal_rng(mu, sigma);
        }
}

//generated quantities {
//  real t_ppc[N_ppc];
//  
//  for (i in 1:N_ppc) {
//    t_ppc[i] = exponential_rng(beta_);
//  }
//}

//data {
//  real N;
//  real sigma;
//}
//
//
//parameters {
//  real x;
//}
//
//
//model {
//  x ~ normal(mu, sigma);
//}