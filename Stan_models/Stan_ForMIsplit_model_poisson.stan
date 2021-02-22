data {
  int<lower=0> J; 
  int<lower=0> N;
  int<lower=1,upper=J> Exploratory[N];
  vector[N] Inonat_z; 
  vector[N] Idwcut_z; 
  vector[N] Iharv_z; 
  int<lower=0> y[N]; 
}
parameters {
  vector[J] a;
  vector[3] b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b_non;
  real<lower=0,upper=100> sigma_b_cut;
  real<lower=0,upper=100> sigma_b_harv;
}
transformed parameters {

   real lp[N];
   
   real <lower=0> mu[N];

   for (i in 1:N) {
     // Linear predictor
     lp[i] = a[Exploratory[i]] + Inonat_z[i] * b[1] + Idwcut_z[i] * b[2] + Iharv_z[i] * b[3];

     // Mean
     mu[i] = exp(lp[i]);
   }
}
model {
  // priors
  sigma_a ~ cauchy(0, 1);
  a ~ normal(mu_a, sigma_a);

  sigma_b_non ~ cauchy(0, 2.5);
  b[1] ~ normal(0, sigma_b_non);
    
  sigma_b_cut ~ cauchy(0, 2.5);
  b[2] ~ normal(0, sigma_b_cut);
    
  sigma_b_harv ~ cauchy(0, 2.5);
  b[3] ~ normal(0, sigma_b_harv);
  
  // likelihood
  y ~ poisson(mu);

}
