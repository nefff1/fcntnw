data {
  int<lower=0> J; 
  int<lower=0> N;
  int<lower=1,upper=J> Exploratory[N];
  vector[N] LUI_z; 
  int<lower=0> y[N]; 
}
parameters {
  vector[J] a;
  real b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b;
}
transformed parameters {

   real lp[N];
   
   real <lower=0> mu[N];

   for (i in 1:N) {
     // Linear predictor
     lp[i] = a[Exploratory[i]] + LUI_z[i] * b;

     // Mean
     mu[i] = exp(lp[i]);
   }
}
model {
  // priors
  sigma_a ~ cauchy(0, 1);
  a ~ normal(mu_a, sigma_a);

  sigma_b ~ cauchy(0, 2.5);
  b ~ normal(0, sigma_b);
  
  // likelihood
  y ~ poisson(mu);

}
