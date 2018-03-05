# bayesian_ICC
R Class to estimate intra class correlations (ICCs) in a Bayesian framework.

Currently only ICC1, ICC1k, ICC3, and ICC3k (ICC2/2k are easy to add, will be done soon).

I use Shrout and Fleiss' (1979) terminology, which is also followed in the psych package.

This class estimate ICC values for many variables at the same time and in the same model; just give it a character vector with all desired variable as the first argument.

## On the test value of 0.01

The hypotheses test, by default, if the estimated ICC value is greater than 0.01. Why not 0, you may ask?

The Bayesian estimation ensures that all estimated ICC values will be in the interval [0,1]. Thus, the confidence intervals for the estimated ICC values will always exclude 0, and this is not a meaningful test. To check this, just run fit$icc1(test_value = 0) on a fit object. All ICCs will be significant, no matter how small.

The value of 0.01 is just a random small number. Feel free to set your own lower bound using the test_value parameter to all icc functions.

## Example using data from the multilevel package

```r
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
```
