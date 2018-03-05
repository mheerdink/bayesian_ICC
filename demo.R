## This script demonstrates all six ICCs and compares them to the ICCs estimated
## using a frequentist approach (in the psych package. 

library(brms)
library(rstan)
library(dplyr)
library(tidyr)
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

# set random seed for reproducibility
# (the icc estimates will still vary from run to run)
set.seed(123)

# get the data in long format for the bayesian ICC and create a
# second rating column to demonstrate the multivariate ICC
sf.df <- sf %>%
    as_data_frame() %>%
    mutate(object = 1:n()) %>%
    gather(judge, rating, -object) %>%
    mutate(rating2 = rnorm(nrow(.)))

head(sf.df, 8)

#######################
### ICC1 and ICC 1k ###
#######################

# create the bICC object with type = 1
b <- bayesian_ICC(c('rating', 'rating2'), 'object', type = 1)

# Fit the model
# You might want to adjust the priors; all priors are uniform on [-∞, ∞] by default, run b$get_priors(df.df)
# to get the priors and later pass them to fit using the 'prior' argument, e.g. b$fit(sf.df, prior = b$get_priors(df.df))
b$fit(sf.df)

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

# now calculate and plot ICC1 and ICC1k
b$icc()
plot(b$icc())

# are the values different from the ICC1 and ICC1k values found by psych::ICC()?
# (look in rows 1 and 3)
b$icc(test = paste('=', c(.17, 0, .44, 0)))
# conclusion: slightly higher

# and, just for the fun of it, something that you cannot do with 'normal' ICC values:
# how much higher are the ICC values for 'rating' than those for the randomly generated 'rating2'?
hypothesis(b$get_fit(), c(paste(b$icc1_formulae()[1:2], collapse=' > '), paste(b$icc1_formulae()[3:4], collapse=' > ')), class=NULL)

#######################
### ICC2 and ICC 2k ###
#######################

# create the bICC object with type = 2
b <- bayesian_ICC(c('rating', 'rating2'), 'object', 'judge', type = 2)

# Fit the model
# (you might want to adjust the default uniform priors, run b$get_priors(df.df) to get a starting point)
b$fit(sf.df, thin = 8, iter = 8000, warmup = 200, control=list(adapt_delta = .99))

# check sampling
stan_trace(b$get_fit()$fit)
stan_ac(b$get_fit()$fit)
# looks good

# now calculate and plot ICC2 and ICC2k
b$icc()
plot(b$icc())

# are the values different from the ICC2 and ICC2k values found by psych::ICC()?
# (look in rows 1 and 3)
b$icc(test = paste('=', c(.29, 0, .62, 0)))
# conclusion: ICC2k is quite a bit lower

#######################
### ICC3 and ICC 3k ###
#######################

# create the bICC object with type = 3
b <- bayesian_ICC(c('rating', 'rating2'), 'object', 'judge', type = 3)

# Fit the model
# (you might want to adjust the default uniform priors; run b$get_priors(df.df) to get the starting point)
b$fit(sf.df, thin = 8, iter = 8000, warmup = 200, control=list(adapt_delta = .99))

# check sampling:
stan_trace(b$get_fit()$fit)
stan_ac(b$get_fit()$fit)
# all good

# now calculate and plot ICC3 and ICC3k
b$icc()
plot(b$icc())

# are the values different from the ICC3 and ICC3k values found by psych::ICC()?
# (look in rows 1 and 3)
b$icc(test = paste('=', c(.71, 0, .91, 0)))
# conclusion: very similar
