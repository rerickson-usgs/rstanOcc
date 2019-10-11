context("Test logistic target function in Stan.")
library(rstanOcc)
options(mc.cores = parallel::detectCores())


set.seed(12345)
n_reps = 10000
y <- rbinom(n = n_reps * 2, 1, c(0.25, 0.75))
x <- rep(c(0, 1), times = n_reps)

l_out <- logistic_target_stan(y, x)

alpha <- rstan::summary(l_out)$summary[ ,"mean"][1]
beta <- alpha + rstan::summary(l_out)$summary[ ,"mean"][2]


test_that("logistic_target_stan works as expected",
{
    expect_equal( as.numeric(c(plogis(alpha), plogis(beta))),
                 aggregate(y, by = list(x), mean)$x, tolerance = 0.001)
}
)
