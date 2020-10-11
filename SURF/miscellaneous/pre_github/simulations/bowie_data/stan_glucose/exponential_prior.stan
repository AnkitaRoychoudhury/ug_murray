data {
  int N;
  real beta_mu;
  real beta_sigma;
}


generated quantities {
  // Parameters
  real<lower=0> beta_;

  // Data
  real t[N];

  beta_ = lognormal_rng(beta_mu, beta_sigma);

  for (i in 1:N) {
    t[i] = exponential_rng(beta_);
  }
}