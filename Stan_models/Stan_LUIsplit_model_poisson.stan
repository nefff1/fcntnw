data {
  int<lower=0> J; 
  int<lower=0> N;
  int<lower=1,upper=J> Exploratory[N];
  vector[N] M_z; 
  vector[N] F_z; 
  vector[N] G_z; 
  int<lower=0> y[N]; 
}
parameters {
  vector[J] a;
  vector[3] b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b_m;
  real<lower=0,upper=100> sigma_b_f;
  real<lower=0,upper=100> sigma_b_g;
}
transformed parameters {

   real lp[N];
   
   real <lower=0> mu[N];

   for (i in 1:N) {
     // Linear predictor
     lp[i] = a[Exploratory[i]] + M_z[i] * b[1] + F_z[i] * b[2] + G_z[i] * b[3];

     // Mean
     mu[i] = exp(lp[i]);
   }
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
  
  // likelihood
  y ~ poisson(mu);

}
