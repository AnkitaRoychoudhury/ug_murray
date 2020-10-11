data {
  int N;
  int N_ppc;
  real t[N];
}

parameters {
  real<lower=0> beta_;
}

model {
  beta_ ~ lognormal(0, 100000);

  t ~ exponential(beta_);
}

generated quantities {
  real t_ppc[N_ppc];
  
  for (i in 1:N_ppc) {
    t_ppc[i] = exponential_rng(beta_);
  }
}