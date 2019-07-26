data {
  // site-level occupancy covariates
  int<lower = 1> nSites;
  int<lower = 1> nPsiCoef;
  matrix[nSites, nPsiCoef] Xpsi;
  
  // sample-level detection covariates
  int<lower = 1> totalSamples;
  int<lower = 1> nPCoef;
  int<lower = 1> nThetaCoef;
  matrix[totalSamples, nPCoef] Vp;
  matrix[totalSamples, nThetaCoef] Wtheta;

  // sample level information  
  int<lower = 0> y[totalSamples];
  int<lower = 0, upper = 1> aObs[totalSamples];
  int<lower = 0> k[totalSamples];
  int<lower = 0, upper = totalSamples> startIndex[nSites];
  int<lower = 0, upper = totalSamples> endIndex[nSites];
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> zObs[nSites];
  
  // number of samples at each site
  int<lower = 0> nSamples[nSites];
}
parameters {
  vector[nPsiCoef] beta_psi;
  vector[nPCoef]   delta_p;
  vector[nThetaCoef]   alpha_theta;
}
transformed parameters {
  vector[totalSamples] logit_p     = Vp     * delta_p;
  vector[totalSamples] logit_theta = Wtheta * alpha_theta;
  vector[nSites] logit_psi         = Xpsi   * beta_psi;
}
model {
  real targetKnownDetection;
  real targetMissedDetectionSample;
  
  vector[nSites] log_psi   = log_inv_logit(logit_psi);
  vector[nSites] log1m_psi = log1m_inv_logit(logit_psi);

  vector[totalSamples] log_theta   = log_inv_logit(logit_theta);
  vector[totalSamples] log1m_theta = log1m_inv_logit(logit_theta);
  
  beta_psi    ~ normal(0, 1);
  alpha_theta ~ normal(0, 1);
  delta_p      ~ normal(0, 1);


  for (site in 1:nSites) { 
    if (nSamples[site] > 0) { // Make sure site was samples
      if (zObs[site] > 0 ) { // Site has known detections
	targetKnownDetection = 0;	
	for(sample in 1:nSamples[site]){ 
	  if(aObs[startIndex[site]:endIndex[site]][sample] > 0){ // Sample has known  detection
	    targetKnownDetection +=
	      log_theta[startIndex[site]:endIndex[site]][sample] +
	      binomial_logit_lpmf(y[startIndex[site]:endIndex[site]][sample] |
				  k[startIndex[site]:endIndex[site]][sample],
				  logit_p[startIndex[site]:endIndex[site]][sample]);
	  } else { // Sample does not have a known detection 
	    targetKnownDetection = log_sum_exp(
				     log_theta[startIndex[site]:endIndex[site]][sample] +
				     binomial_logit_lpmf(y[startIndex[site]:endIndex[site]][sample] |
							 k[startIndex[site]:endIndex[site]][sample],
							 logit_p[startIndex[site]:endIndex[site]][sample]),
				     log1m_theta[startIndex[site]:endIndex[site]][sample]);
	  }
	}
	target += log_psi[site] + targetKnownDetection;
      } else {
	targetMissedDetectionSample = 0;
        // site may or may not be occupied 
	for(sample in 1:nSamples[site]){
	  targetMissedDetectionSample += // Missed with detection or sample
	    log_sum_exp(
			log_theta[startIndex[site]:endIndex[site]][sample] +
			binomial_logit_lpmf(y[startIndex[site]:endIndex[site]][sample] |
					    k[startIndex[site]:endIndex[site]][sample],
					    logit_p[startIndex[site]:endIndex[site]][sample]),
			log1m_theta[startIndex[site]:endIndex[site]][sample]);
	}
	// log1m_psi is probability site is not occupied
	target += log_sum_exp( log_psi[site] + targetMissedDetectionSample,  log1m_psi[site]);
      }
    }
  }  

}

