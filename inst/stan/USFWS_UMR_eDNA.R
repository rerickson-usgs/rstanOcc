data {
  // number of total observation
  int<lower = 1> n_obs;
  
  // site-level covariates
  int<lower = 1> n_sites;
  int<lower = 1> n_psi_coef;
  matrix[ n_sites, n_psi_coef] X_psi;
  int<lower = 0> site_detection[n_sites];

  
  // sample-level capture covariates
  int<lower = 1> n_samples[n_sites];
  int<lower = 1> n_theta_coef;
  matrix[ n_obs, n_theta_coef] W_theta;

  // molecular-level detection covariates
  int<lower = 1> n_p_coef;
  matrix[ n_obs, n_p_coef] V_p;
  int<lower = 0> AC1[n_obs];
  int<lower = 0> AC3[n_obs];
  
  // sample-level information
  int<lower = 0, upper = n_obs> start_index[n_sites];
  int<lower = 0, upper = n_obs> end_index[n_sites];

  // summary if any samples tested posistive
  int<lower = 0, upper = 1> any_detected[n_obs];

  // molecular replicates per sample 
  int<lower = 1> k[n_obs];

  // create summary input matrix
  int<lower = 1> n_row_sep;
  int<lower = 1> n_col_sep;
  matrix[ n_row_sep, n_col_sep] site_event_predict;

  int<lower = 1> n_predict_samples;

  // priors
  real delta_sd;
  real alpha_sd;
  real beta_sd;

  real delta_mean;
  real alpha_mean;
  real beta_mean;
}
parameters {
  vector[n_p_coef]     delta_p_AC1;
  vector[n_p_coef]     delta_p_AC3;
  vector[n_theta_coef] alpha_theta;
  vector[n_psi_coef]   beta_psi;
}
transformed parameters {
  
  vector[n_obs] logit_p_AC1 = V_p     * delta_p_AC1;
  vector[n_obs] logit_p_AC3 = V_p     * delta_p_AC3;
  vector[n_obs] logit_theta = W_theta * alpha_theta;
  vector[n_sites] logit_psi = X_psi   * beta_psi;
}
model {
  vector[n_obs] log_theta   = log_inv_logit(logit_theta);
  vector[n_obs] log1m_theta = log1m_inv_logit(logit_theta);

  vector[n_sites] log_psi   = log_inv_logit(logit_psi);
  vector[n_sites] log1m_psi = log1m_inv_logit(logit_psi);

  real target_temp;
  
  delta_p_AC1 ~ normal(delta_mean, delta_sd);
  delta_p_AC3 ~ normal(delta_mean, delta_sd);
  alpha_theta ~ normal(alpha_mean, alpha_sd);
  beta_psi    ~ normal(beta_mean, beta_sd);

  target_temp = 0;
 
  for (site in 1:n_sites){
    if(site_detection[site] > 0){
      // Yes, DNA found at site
      for (sample in 1:n_samples[site]){
	  if ( any_detected[start_index[site]:end_index[site]][sample] > 0){
	    // DNA detected in Sample
	    target_temp +=  log_theta[start_index[site]:end_index[site]][sample] +
	      // prob DNA in sample
	      binomial_logit_lpmf( AC1[start_index[site]:end_index[site]][sample] |
				   k[start_index[site]:end_index[site]][sample],
				   logit_p_AC1[start_index[site]:end_index[site]][sample]) +
	      binomial_logit_lpmf( AC3[start_index[site]:end_index[site]][sample] |
				   k[start_index[site]:end_index[site]][sample],
				   logit_p_AC3[start_index[site]:end_index[site]][sample]);
	    // prob of DNA detected in a given sample
	  } else {
	    // no detections in sample
	    target_temp +=
	      log_sum_exp(
			  log_theta[start_index[site]:end_index[site]][sample] +
			  // prob in sample, but ....
			  binomial_logit_lpmf( AC1[start_index[site]:end_index[site]][sample] |
					       k[start_index[site]:end_index[site]][sample],
					       logit_p_AC1[start_index[site]:end_index[site]][sample] ) +
			  binomial_logit_lpmf( AC3[start_index[site]:end_index[site]][sample] |
					       k[start_index[site]:end_index[site]][sample],
					       logit_p_AC3[start_index[site]:end_index[site]][sample] ),
			  // missed DNA with molecular method 
			  log1m_theta[start_index[site]:end_index[site]][sample]
			  // or DNA was simply not in sample
			  );
	  }	
      }
      // end of site known. Add probs togehter and reset temp targets 
      target += target_temp + logit_psi[site]; 
      target_temp = 0;
    } else {
      // no DNA observed at site
      for (sample in 1:n_samples[site]){
	// go through samples 
	target_temp +=
	  log_sum_exp(
		      log_theta[start_index[site]:end_index[site]][sample] +
		      // prob in sample, but ....
		      binomial_logit_lpmf( AC1[start_index[site]:end_index[site]][sample] |
					   k[start_index[site]:end_index[site]][sample],
					   logit_p_AC1[start_index[site]:end_index[site]][sample] ) +
		      binomial_logit_lpmf( AC3[start_index[site]:end_index[site]][sample] |
					   k[start_index[site]:end_index[site]][sample],
					   logit_p_AC3[start_index[site]:end_index[site]][sample] ),
		      // missed DNA with molecular method 
		      log1m_theta[start_index[site]:end_index[site]][sample]
		      // or DNA was simply not in sample
		      );
      }
      // add target probs here
      target += log_sum_exp(log1m_psi[site], logit_psi[site] +  target_temp);
      target_temp = 0;
    }
  } // end of site loop 
}
generated quantities {
  vector<lower = 0, upper = 1>[n_psi_coef] psi_site;
  vector<lower = 0, upper = 1>[n_sites] theta_site;
  vector<lower = 0, upper = 1>[n_sites] p_AC1_site; 
  vector<lower = 0, upper = 1>[n_sites] p_AC3_site; 
  vector<lower = 0, upper = 1>[n_sites] p_posistive;
  
  matrix[n_sites, n_predict_samples] prob_detect_one;

  for(site in 1:n_psi_coef){
    psi_site[site] = inv_logit(beta_psi[site]);
  }

  theta_site = inv_logit(site_event_predict * alpha_theta);
  p_AC1_site = inv_logit(site_event_predict * delta_p_AC1);
  p_AC3_site = inv_logit(site_event_predict * delta_p_AC3);

  for(site in 1:n_sites){
    if(site_detection[site] > 0){
      p_posistive[site] =
  	theta_site[site] *
  	(1.0 - ( 1.0 - p_AC1_site[site]) ^ 8.0) *
  	(1.0 - ( 1.0 - p_AC3_site[site]) ^ 8.0);
    } else {
      p_posistive[site] = 0;
  	}

    for( rep in 1:n_predict_samples){
      prob_detect_one[site, rep] = 1.0 - (1.0 - p_posistive[site]) ^ rep;
    }
  }
}
