#' Bayesian logistic regression using traget probabiltiy with Stan
#'
#' This differes from \code{logistic_stan()} becaues of how the
#' Stan code is written. Both this function and that function
#' are included as teaching tools to build up towards occupancy models.
#' 
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' 2+2 
#' @export
logistic_target_stan <- function(x, y, ...) {
  standata <- list(x = x, y = y, N = length(y))
  out <- rstan::sampling(stanmodels$logistic__target_regression,
                         data = standata, ...)
  return(out)
}
