
data{
  int n;
  int nobs;
  real xobs[nobs];
  real y[n];
  int indx[nobs];
}

parameters {
  real alpha;
  real beta;
  real sigmay;
  real sigmax;
  real x[n];
}

model {
  // priors
  alpha ~ normal(0,100);
  beta ~ normal(0,100);
  sigmay ~ uniform(0,1000);
  sigmax ~ uniform(0,1000);
  
  // model structure  
  for (i in 1:nobs){
    xobs[i] ~ normal(x[indx[i]], sigmax);
  }
  for (i in 1:n){
    y[i] ~ normal(alpha + beta*x[i], sigmay);
  }
}
  