data {
  int<lower=0> J; 
  int<lower=0> N;
  int<lower=1,upper=J> Exploratory[N];
  vector[N] M_z; 
  vector[N] F_z; 
  vector[N] G_z; 
  vector[N] y; 
}
parameters {
  vector[J] a;
  vector[3] b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b_m;
  real<lower=0,upper=100> sigma_b_f;
  real<lower=0,upper=100> sigma_b_g;
  real<lower=0,upper=100> sigma_y;
}
transformed parameters {
    vector[N] y_hat;

  for (i in 1:N)
  y_hat[i] = a[Exploratory[i]] + M_z[i] * b[1] + F_z[i] * b[2] + G_z[i] * b[3];
}
model {
  // priors
  sigma_a ~ cauchy(0, 1);
  a ~ normal(mu_a, sigma_a);

  sigma_b_m ~ cauchy(0, 2.5);
  b[1] ~ normal(0, sigma_b_m);
    
  sigma_b_f ~ cauchy(0, 2.5);
  b[2] ~ normal(0, sigma_b_f);
    
  sigma_b_g ~ cauchy(0, 2.5);
  b[3] ~ normal(0, sigma_b_g);

  sigma_y ~ cauchy(0, 25);
  
  // likelihood
  y ~ normal(y_hat, sigma_y);

}
