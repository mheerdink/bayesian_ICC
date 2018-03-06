# bayesian_ICC

R Class to fit intra class correlations (ICCs) on multiple variables simultaneously using a Bayesian approach.

It can estimate ICC values for many variables at the same time and in the same model; just give it a character vector with all desired variable as the first argument. This allows you to do some fancy things like studying the difference between ICCs for different variables.

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

set.seed(123)

# the data have to be in long format for the bayesian ICC
# also, create a second rating column to demonstrate the multivariate ICC
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
```

## Caveats

Bayesian estimation of the confidence intervals of ICCs is not perfect, as explained in this paper: https://doi.org/10.1186/1471-2288-14-121

The authors recommend the present approach only if the data are normally distributed, and if the number of levels of each random factors (i.e., the number of objects or the number of objects) is greater than 8. In other cases, Modified Large Sample (MLS) or Generalized Confidence Interval (GCI) methods might be better for obtaining the confidence interval. This means that, strictly speaking, this approach cannot be used with the Shrout & Fleiss (1979) example, as I did above.

## About priors

Additionally, the same authors recommend (based on Gelman's work) a uniform flat (improper) prior on [0, ∞]. However, this is against the recommendations on priors on the Stan website (https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations), I find that this often leads to problems in estimating the model, and I can't seem to figure out how to set the uniform prior to be one-sided.

The default priors are therefore cauchy(0, 3) priors, which most likely cover the likely locations of the sd for the intercept well (especially if working with z-scaled predictors).

You can also use either:

* BRMS default priors (currently student-t): call `p <- b$get_priors(df, priors = 'brms')` followed by `b$fit(df, prior=p)`
* Uniform priors on [-∞, ∞]: call `p <- b$get_priors(df, priors = 'uniform')` followed by `b$fit(df, prior=p)`

I hope this will be useful to someone!

## Acknowledgements

None of this is my own invention, of course. All I did was combining others' work into an easy to use package.

I borrowed formulae, inspiration, and programming ideas from:

* http://www.kohleth.info/2017/10/09/all-the-different-iccs/
* https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
* http://www.cyclismo.org/tutorial/R/s3Classes.html
* Murray Logan's excellent course on both Frequentist and Bayesian statistics, which I can highly recommend, http://www.flutterbys.com.au/stats/course.html
* Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in assessing rater reliability. *Psychological Bulletin*, *86*, 420–428. https://doi.org/10.1037/0033-2909.86.2.420

Credits also go to the psych, brms, dplyr/tidyr, and stan creators.
