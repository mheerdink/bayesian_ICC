library(brms)
library(rstan)
library(dplyr)
library(psych)

# Use example from Shrout and Fleiss (1979), also used in the psych::ICC() help page
sf <- matrix(c(
    9L,    2L,   5L,    8L,
    6L,    1L,   3L,    2L,
    8L,    4L,   6L,    8L,
    7L,    1L,   2L,    6L,
    10L,   5L,   6L,    9L,
    6L,    2L,   4L,    7L), ncol=4, byrow=TRUE)
colnames(sf) <- paste("J", 1:4, sep="")
rownames(sf) <- paste("S", 1:6, sep="")
ICC(sf)

source('bayesian_ICC.R')

# get the data in long format for the bayesian ICC
sf.df <- sf %>%
    as_data_frame() %>%
    mutate(object = 1:n()) %>%
    gather(judge, value, -object)

sf.df

#######################
### ICC1 and ICC 1k ###
#######################

b <- bayesian_ICC('value', 'object', which = 1)

# Fit the model (you might want to adjust the priors)
b$fit(sf.df, prior = c(prior(cauchy(0, 5), class='sigma')))

# adapt_delta should be increased
# in addition, it is advisable to check if sampling went allright
# by exploring the sampler statistics, e.g.: (using the rstan package)
stan_trace(b$get_fit()$fit)
stan_ac(b$get_fit()$fit)
# mixing is fine, but there's some autocorrelation

# refit with higher adapt_delta and thinning:
b$refit(thin = 8, iter = 8000, warmup = 200, control = list(adapt_delta = .99))
stan_ac(b$get_fit()$fit)
# no more issues

# get ICC1
b$icc1()
plot(b$icc1())
# is it similar to the ICC1 value found by psych::ICC()?
b$icc1(test_value = .17)
# conclusion: not different

# get ICC1k
b$icc1k()
plot(b$icc1k())
# is it similar to the ICC1k value found by psych::ICC()?
b$icc1k(test_value = .44)
# conclusion: not different

#######################
### ICC3 and ICC 3k ###
#######################

b <- bayesian_ICC('value', 'object', 'judge', which = 3)

# Fit the model (you might want to adjust the priors)
b$fit(sf.df, prior = c(prior(cauchy(0, 5), class='sigma')))

# adapt_delta should be increased
# in addition, it is advisable to check if sampling went allright
# by exploring the sampler statistics, e.g.: (using the rstan package)
stan_trace(b$get_fit()$fit)
stan_ac(b$get_fit()$fit)
# mixing is fine, but there's some autocorrelation

# refit with higher adapt_delta and thinning:
b$refit(thin = 8, iter = 8000, warmup = 200, control = list(adapt_delta = .99))
stan_ac(b$get_fit()$fit)
# no more issues

# get ICC3
b$icc3()
plot(b$icc3())
# is it similar to the ICC3 value found by psych::ICC()?
b$icc3(test_value = .71)
# conclusion: not different

# get ICC1k
b$icc1k()
plot(b$icc1k())
# is it similar to the ICC3k value found by psych::ICC()?
b$icc3k(test_value = .91)
# conclusion: not different