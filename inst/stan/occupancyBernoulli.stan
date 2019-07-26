data {
  // site-level occupancy covariates
  int<lower = 1> nSites;
  int<lower = 1> nPsiCoef;
  matrix[nSites, nPsiCoef] Xpsi;
  
  // survey-level detection covariates
  int<lower = 1> totalSurveys;
  int<lower = 1> nPCoef;
  matrix[totalSurveys, nPCoef] Vp;

  // survey level information  
  int<lower = 1, upper = nSites> site[totalSurveys];
  int<lower = 0, upper = 1> y[totalSurveys];
  int<lower = 0, upper = totalSurveys> startIndex[nSites];
  int<lower = 0, upper = totalSurveys> endIndex[nSites];
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> z[nSites];
  
  // number of surveys at each site
  int<lower = 0> nSurveys[nSites];
}
parameters {
  vector[nPsiCoef] beta_psi;
  vector[nPCoef] beta_p;
}
transformed parameters {
  vector[totalSurveys] logit_p = Vp * beta_p;
  vector[nSites] logit_psi = Xpsi * beta_psi;
}
model {
  vector[nSites] log_psi = log_inv_logit(logit_psi);
  vector[nSites] log1m_psi = log1m_inv_logit(logit_psi);
  
  beta_psi ~ normal(0, 1);
  beta_p ~ normal(0, 1);
  for (i in 1:nSites) {
    if (nSurveys[i] > 0) {
      if (z[i] > 0 ) {
        // site is occupied
        target += log_psi[i] + 
	  bernoulli_logit_lpmf(y[startIndex[i]:endIndex[i]] | 
			       logit_p[startIndex[i]:endIndex[i]]);
      } else {
        // site may or may not be occupied
        target += log_sum_exp(
			      log_psi[i] +
			      bernoulli_logit_lpmf(y[startIndex[i]:endIndex[i]] |
						   logit_p[startIndex[i]:endIndex[i]]), 
			      log1m_psi[i]);
      }
    }
  }
}
