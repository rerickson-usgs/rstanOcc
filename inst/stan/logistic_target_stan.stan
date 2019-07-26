data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real beta;
}
model {
  
  for( index in 1:N){
    target += binomial_logit_lpmf( y[index] | 1, alpha + beta * x[index]); 
      }
}

