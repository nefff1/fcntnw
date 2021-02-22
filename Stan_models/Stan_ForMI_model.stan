data {
  int<lower=0> J; 
  int<lower=0> N;
  int<lower=1,upper=J> Exploratory[N];
  vector[N] ForMI_z; 
  vector[N] y; 
}
parameters {
  vector[J] a;
  real b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b;
  real<lower=0,upper=100> sigma_y;
}
transformed parameters {
    vector[N] y_hat;

  for (i in 1:N)
  y_hat[i] = a[Exploratory[i]] + ForMI_z[i] * b;
}
model {
  // priors
  sigma_a ~ cauchy(0, 1);
  a ~ normal(mu_a, sigma_a);

  sigma_b ~ cauchy(0, 2.5);
  b ~ normal(0, sigma_b);

  sigma_y ~ cauchy(0, 25);
  
  // likelihood
  y ~ normal(y_hat, sigma_y);

}
