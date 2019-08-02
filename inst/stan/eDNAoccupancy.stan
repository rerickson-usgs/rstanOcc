data {
  // site-level occupancy covariates
  int<lower = 1> n_sites;
  int<lower = 1> n_Psi_coef;
  matrix[n_sites, n_Psi_coef] X_psi;
  
  // sample-level detection covariates
  int<lower = 1> n_samples;
  int<lower = 1> n_p_coef;
  int<lower = 1> n_Theta_coef;
  matrix[n_samples, n_p_coef] V_p;
  matrix[n_samples, n_Theta_coef] W_theta;

  // sample level information  
  int<lower = 0> y[n_samples];
  int<lower = 0, upper = 1> a_obs[n_samples];
  int<lower = 0> k[n_samples];
  int<lower = 0, upper = n_samples> start_index[n_sites];
  int<lower = 0, upper = n_samples> end_index[n_sites];
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> z_obs[n_sites];
  
  // number of samples at each site
  int<lower = 0> n_samples_per_site[n_sites];
}
parameters {
  vector[n_Psi_coef] beta_psi;
  vector[n_p_coef]   delta_p;
  vector[n_Theta_coef]   alpha_theta;
}
transformed parameters {
  vector[n_samples] logit_p     = V_p     * delta_p;
  vector[n_samples] logit_theta = W_theta * alpha_theta;
  vector[n_sites] logit_psi         = X_psi   * beta_psi;
}
model {
  real targetKnownDetection;
  real targetMissedDetectionSample;
  
  vector[n_sites] log_psi   = log_inv_logit(logit_psi);
  vector[n_sites] log1m_psi = log1m_inv_logit(logit_psi);

  vector[n_samples] log_theta   = log_inv_logit(logit_theta);
  vector[n_samples] log1m_theta = log1m_inv_logit(logit_theta);
  
  beta_psi    ~ normal(0, 1);
  alpha_theta ~ normal(0, 1);
  delta_p      ~ normal(0, 1);


  for (site in 1:n_sites) { 
    if (n_samples_per_site[site] > 0) { // Make sure site was samples
      if (z_obs[site] > 0 ) { // Site has known detections
	targetKnownDetection = 0;	
	for(sample in 1:n_samples_per_site[site]){ 
	  if(a_obs[start_index[site]:end_index[site]][sample] > 0){ // Sample has known  detection
	    targetKnownDetection +=
	      log_theta[start_index[site]:end_index[site]][sample] +
	      binomial_logit_lpmf(y[start_index[site]:end_index[site]][sample] |
				  k[start_index[site]:end_index[site]][sample],
				  logit_p[start_index[site]:end_index[site]][sample]);
	  } else { // Sample does not have a known detection 
	    targetKnownDetection = log_sum_exp(
				     log_theta[start_index[site]:end_index[site]][sample] +
				     binomial_logit_lpmf(y[start_index[site]:end_index[site]][sample] |
							 k[start_index[site]:end_index[site]][sample],
							 logit_p[start_index[site]:end_index[site]][sample]),
				     log1m_theta[start_index[site]:end_index[site]][sample]);
	  }
	}
	target += log_psi[site] + targetKnownDetection;
      } else {
	targetMissedDetectionSample = 0;
        // site may or may not be occupied 
	for(sample in 1:n_samples_per_site[site]){
	  targetMissedDetectionSample += // Missed with detection or sample
	    log_sum_exp(
			log_theta[start_index[site]:end_index[site]][sample] +
			binomial_logit_lpmf(y[start_index[site]:end_index[site]][sample] |
					    k[start_index[site]:end_index[site]][sample],
					    logit_p[start_index[site]:end_index[site]][sample]),
			log1m_theta[start_index[site]:end_index[site]][sample]);
	}
	// log1m_psi is probability site is not occupied
	target += log_sum_exp( log_psi[site] + targetMissedDetectionSample,  log1m_psi[site]);
      }
    }
  }  

}

