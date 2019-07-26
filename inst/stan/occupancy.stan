/**
 * Occupancy model; Single season, no covariates
 *
 * Translated from Marc Kery's WinBUGS model:
 * http://www.fisheriesstockassessment.com/TikiWiki/tiki-index.php?page=BUGS+Occupancy
 *
 * Based on: 
 * J. Andrew Royle and Marc Kery. 2007. A Bayesian state-space
 * formulation of dynamic occupancy models.  Ecology 88(7),
 * 1813--1823.
 * URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3892356/
 */

data {
  int<lower=0> n_sampling_units;
  int<lower=0> n_surveys;
  int<lower=0,upper=1> y[n_sampling_units, n_surveys];
}
parameters {
  real<lower=0,upper=1> psi;
  real<lower=0,upper=1> p;
}
model {
  // local variables to avoid recomputing log(psi) and log(1 - psi)
  real log_psi;
  real log1m_psi;
  log_psi = log(psi);
  log1m_psi = log1m(psi);

  // priors
  psi ~ uniform(0,1);
  p   ~ uniform(0,1);
  
  // likelihood
  for (r in 1:n_sampling_units) {
    if (sum(y[r]) > 0)
      target += log_psi + bernoulli_lpmf(y[r] | p);
    else
      target += log_sum_exp(log_psi + bernoulli_lpmf(y[r] | p),
			    log1m_psi);;
  }
}
