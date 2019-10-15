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
#' @param W_theta matrix of theta predictors
#' @param y vector of number of detections per sample
#' @param a_obs vetor if detections occured within sample MTI
#' @param k vector of number of molecular replicates per sample
#' @param start_index starting index for each site used MTI
#' @param end_index ending index for each site used MTI
#' @param z_obs vector for each site if any detections occured
#' @param n_samples_per_site number of samples per each site 
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
options(mc.cores = parallel::detectCores())


library(dplyr)
set.seed(123456)
## Simulate sites 
n_sites <- 20
psi <- 0.75
z_known <- rbinom(n_sites, size = 1, prob = psi)
site_data <- data.frame(site_id = 1:n_sites,
                        z_known)
X_psi <- model.matrix( ~ 1, data = site_data)
n_Psi_coef <- ncol(X_psi)

## Simulate sampels and sub-samples at each site
n_samples_per_site <- 20
theta <- 0.75
k_all <- 8
p <- 0.75

sample_data <- data.frame(site_id = rep(1:n_sites, each = n_sites),
                          k = k_all)
sample_data <- sample_data %>% full_join(site_data)
sample_data <- mutate(sample_data, theta =theta, p = 0.75)
sample_data <- mutate(sample_data,
                      a_known = z_known * rbinom(nrow(sample_data), 1, theta))
sample_data <- mutate(sample_data,
                      y_known = a_known * rbinom(nrow(sample_data), k, p))
sample_data <- mutate(sample_data,
                      a_obs = ifelse( y_known > 0, 1, 0))
sample_data <-
    sample_data %>%
    mutate(sample_id = 1:nrow(sample_data))
    

## for some scenarios a_obs may not be equal to a_known
## sample_data %>% filter(a_known != a_obs)

## Create predictor matricies
V_p <- model.matrix( ~ 1, sample_data)
n_p_coef <- ncol(V_p)

sample_data <- mutate(sample_data, site_id = factor(site_id))
W_theta <- model.matrix( ~ 1, sample_data)
n_Theta_coef <- ncol(W_theta)

## Calculate "observed" values for checking model
## and our intuition of the model 
site_level_summary <-
    sample_data %>%
    group_by(site_id) %>%
    summarize(
        theta_obs = mean(a_known),
        n_samples_per_site = n()
    ) %>%
    mutate(z_obs = ifelse(theta_obs > 0, 1, 0))

theta_obs <-
    sample_data %>%
    filter(z_known > 0) %>%
    pull(a_known) %>%
    mean(.)

psi_obs <- site_level_summary %>% pull(z_obs) %>% mean(.)
p_obs <-
    sample_data %>%
    filter(z_known > 0 & a_known > 0) %>%
    pull(y_known) %>%
    mean(.) / 8

## Create sample index
start_end_index <-
    sample_data %>%
    group_by(site_id) %>%
    summarize(start_index = min(sample_id),
              end_index = max(sample_id))

## Create inputs for model 
eDNA_model_out <-
    eDNAoccupancy_stan(
    n_sites = n_sites,
    n_Psi_coef = n_Psi_coef,
    X_psi = X_psi,
    n_samples = nrow(sample_data),
    n_p_coef = n_p_coef,
    n_Theta_coef = n_Theta_coef,
    V_p = V_p,
    W_theta = W_theta,
    y = sample_data %>% pull(y_known),
    a_obs = sample_data %>% pull(a_obs),
    k = sample_data %>% pull(k),
    start_index = start_end_index %>% pull(start_index),
    end_index = start_end_index %>% pull(end_index),
    z_obs = site_level_summary %>% pull(z_obs),
    n_samples_per_site = site_level_summary %>% pull(n_samples_per_site))

p_obs
theta_obs
psi_obs


compare <-
    rbind(c(psi_obs, theta_obs, p_obs), 
          plogis(rstan::summary(eDNA_model_out,
                                pars = c("beta_psi", "alpha_theta", "delta_p"))$summary[ , 'mean']))

#' @export
eDNAoccupancy_stan <- 
    function(n_sites, n_Psi_coef, X_psi,
             n_samples, n_p_coef, n_Theta_coef,
             V_p, W_theta, y, a_obs, k,
             start_index, end_index,
             z_obs, n_samples_per_site,  ...){

        stan_data <- list(
            n_sites, n_Psi_coef, X_psi,
            n_samples, n_p_coef, n_Theta_coef,
            V_p, W_theta, y, a_obs, k,
            start_index, end_index,
            z_obs, n_samples_per_site)

        out <- rstan::sampling(stanmodels$eDNAoccupancy,
                               data = stan_data, ...)
        return(out)
    }
