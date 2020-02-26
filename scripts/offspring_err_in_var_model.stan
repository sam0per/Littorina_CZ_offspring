//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data are vectors 'x_meas', 'tau' and 'y_mean' of length 'N'.
data {
  int<lower=0> N;
  vector[N] x_meas;
  real<lower=0> tau[N];
  vector[N] y_mean;
}

// The parameters accepted by the model.
parameters {
  real x[N];
  real mu_x;
  real sigma_x;
  real alpha;
  real beta;
  real<lower=0> sigma_y;
}

// The model to be estimated. We model the output
// 'y_mean' to be normally distributed with mean 'alpha + beta * x'
// and standard deviation 'sigma_y'.
model {
  x ~ normal(mu_x, sigma_x);
  x_meas ~ normal(x, tau);
  alpha ~ normal(0, 5);
  beta ~ normal(1, 5);
  sigma_y ~ cauchy(0, 5);
  y_mean ~ normal(alpha + beta * x, sigma_y);
}

