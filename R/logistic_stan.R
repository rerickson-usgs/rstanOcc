#' Simple bayesian logistic regression with Stan
#'
#' 
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' set.seed(12345)
#' n_reps = 10000
#' y <- rbinom(n = n_reps * 2, 1, c(0.25, 0.75))
#' x <- rep(c(0, 1), times = n_reps)
#'
#' l_out <- logistic_stan(y, x)
#'
#' ## show simulated data is almost the same as the recovered probs
#' print(rstan::summary(l_out)$summary[ ,"mean"][1:2])
#' print(aggregate(y, by = list(x), mean)$x)
#' 
#' 

#' @export
logistic_stan <- function(x, y, ...) {
    standata <- list(x = x, y = y, N = length(y))
    out <- rstan::sampling(stanmodels$logistic_stan, data = standata, ...)
    return(out)
}
