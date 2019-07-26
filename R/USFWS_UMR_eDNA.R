#' Occupancy model used by USFWS for eDNA on UMR
#'
#' This model is was used by the USFWS on the UMR as part
#' of their bighead and silver carp monitoring.
#'
#' This will include a reference to Mize et al. once the
#' help file is fleshed out. 
#'
#' With this model, I will need to move more details into hte r
#' 
#' @param n_obs number of total observations
#' @param n_sites number of sampling sites (e.g., pools)
#' @param n_psi_coef number of site-level coefficients
#' @param X_psi matrix for site-level predictors
#' @param site_detections binar vector of eDNA detetion at each site
#' @param n_samples vector of no. samples per site
#' @param n_theta_coef number of sample-level coeffients
#' @param W_theta matrix for sample-level predictors
#' @param n_p_coef number of detection-level coeffients
#' @param V_p matrix for sample-level predictors
#' @param AC1 number of positive detections per sample using ACTM 1 assasy
#' @param AC3 number of positive detections per sample using ACTM 3 assasy
#' @param start_index sample-level information
#' @param end_index sample-level information
#' @param any_detected were there any_detected in sample? 0 or 1
#' @param k number of molecular replicates per sample
#' @param n_row_sep number of rows in predictor matrix
#' @param n_col_sep number of columns in predicotr matrix
#' @param site_event_predict matrix of site-event predictions
#' @param n_predict_samples number of prediciton samples
#' @param delta_sd prior for delta SD
#' @param alpha_sd prior for alpha SD
#' @param beta_sd prior for beta SD
#' @param delta_mean prior for delta mean
#' @param alpha_mean prior for alpha mean
#' @param beta_mean prior for beta mean 
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' 2 + 2
#' @export
USFWS_UMR_eDNA_stan <- function(n_obs, n_sites, n_psi_coef, X_psi,
                                site_detections, n_samples, n_theta_coef,
                                V_p, AC1, AC3, start_index, end_index,
                                any_detected, k, n_row_sep, n_col_sep,
                                site_event_predict, n_predict_samples,
                                delta_sd, alpha_sd, beta_sd,
                                delta_mean, alpha_mean, beta_mean, ...) {
  standata <- list(n_obs, n_sites, n_psi_coef, X_psi,
                                site_detections, n_samples, n_theta_coef,
                                V_p, AC1, AC3, start_index, end_index,
                                any_detected, k, n_row_sep, n_col_sep,
                                site_event_predict, n_predict_samples,
                                delta_sd, alpha_sd, beta_sd,
                                delta_mean, alpha_mean, beta_mean)
  out <- rstan::sampling(stanmodels$logistic_regression, data = standata, ...)
  return(out)
}
