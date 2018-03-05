# bayesian_ICC

R Class to fit intra class correlations (ICCs) on multiple variables simultaneously using a Bayesian approach.

It can estimate ICC values for many variables at the same time and in the same model; just give it a character vector with all desired variable as the first argument. This allows you to do some fancy things like comparing if the ICC of one variable is different from the same ICC for a different variable.

In Shrout and Fleiss' (1979) terminology (which is also followed in the psych package), it estimates ICC1, ICC1k, ICC2, ICC2k, ICC3, and ICC3k.

## Example usage

```r
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

# get the data in long format for the bayesian ICC
# and create a second rating column to demonstrate the multivariate ICC
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
# (you might want to adjust the priors; to see the default priors, run b$get_priors(df.df))
b$fit(sf.df, prior = c(prior(cauchy(0, 5), class='sigma', resp='rating'), prior(cauchy(0, 5), class='sigma', resp='rating2')))

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
# conclusion: quite similar

# and, just for the fun of it, something that you cannot do with 'normal' ICC values:
# how much higher are the ICC values for 'rating' than those for the randomly generated 'rating2'?
hypothesis(b$get_fit(), c(paste(b$icc1_formulae()[1:2], collapse=' > '), paste(b$icc1_formulae()[3:4], collapse=' > ')), class=NULL)
```

## Acknowledgements

None of this is my own invention, of course. All I did was combining others' work into an easy to use package.

I borrowed formulae, inspiration, and programming ideas from:

* http://www.kohleth.info/2017/10/09/all-the-different-iccs/
* https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
* http://www.cyclismo.org/tutorial/R/s3Classes.html
* Murray Logan's excellent course on both Frequentist and Bayesian statistics, which I can highly recommend, http://www.flutterbys.com.au/stats/course.html
* Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in assessing rater reliability. *Psychological Bulletin*, *86*, 420–428. https://doi.org/10.1037/0033-2909.86.2.420

Credits also go to the psych, brms, dplyr/tidyr, and stan creators.
