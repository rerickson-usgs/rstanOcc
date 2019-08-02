#' Three-level occupancy model
#'
#' Parameter with MTI will likely be moved to internal functions
#' in future version of this function 
#'
#' Also, site needs to be changed to sampling event
#' 
#' @param n_sites Number of sites
#' @param n_Psi_coef Number of Psi coefficients
#' @param X_psi matrix of psi predictors
#' @param n_samples total number of samples
#' @param n_p_coef number of p coefficients
#' @param n_Theta_coef number of theta coefficients
#' @param V_p matrix of p predictors
#' @param W_p matrix of theta predictors
#' @param y vector of number of detections per sample
#' @param a_obs vetor if detections occured within sample MTI
#' @param k number of molecular replicates per sample
#' @param start_index starting index for each site used MTI
#' @param end_index ending index for each site used MTI
#' @param z_obs vector for each site if any detections occured
#' @param n_samples_per_site number of samples per each site 
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`

#' @export
eDNAoccupancy_stan <- 
    function(n_sites, n_Psi_coef, X_psi,
             n_samples, n_p_coef, n_Theta_coef,
             V_p, W_p, y, a_obs, k,
             start_index, end_index,
             z_obs,  ...){

        stan_data <- list(
            n_sites, n_Psi_coef, X_psi,
            n_samples, n_p_coef, n_Theta_coef,
            V_p, W_p, y, a_obs, k,
            start_index, end_index,
            z_obs)

        out <- rstan::sampling(stanmodels$eDNAoccupancy,
                               data = stan_data, ...)
        return(out)
    }
