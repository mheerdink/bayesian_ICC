bayesian_ICC <- function(vars, objects, judges = NULL, type = 1, rescor = F) {
    # Class to fit six ICCs (ICC1, ICC1k, ICC2, ICC2k, ICC3, ICC3k) on multiple variables simultaneously in a Bayesian framework with ease
    #
    # Args:
    #  vars      vector of characters, specifying the variable names of the variables to aggregate
    #  objects   character, specifying the objects having been rated (e.g., the number of groups or video clips)
    #  judges    character, optional specification for specification of judges (required to calculate ICC2 and ICC3)
    #  type      ICC type, 1 (each judge judged only one object), 2 (judges judged multiple objects and were randomly sampled) or 3 (judges judged multiple object and are the complete population of possible judges)
    #  rescor    Whether or not to model response correlation if > 1 var
    #
    # Returns:
    #  A class with multiple methods
    
    if (! is.character(vars) | ! is.character(objects) | length(objects) > 1)
        stop("Vars should be a character vector; objects should be a character variable")
    if (! is.null(judges) & (!is.character(judges) | length(judges) > 1))
        stop("Judges should be NULL or a character variable")
    if (length(type) != 1 | ! type %in% 1:3)
        stop("type should be 1, 2, or 3 for ICC1, ICC2, and ICC3, respectively")
    if (type %in% 2:3 & is.null(judges))
        stop("Specify judges for ICC2 or ICC3")
    if (type == 1 & !is.null(judges))
        warning("Judges argument will be ignored when calculating ICC1")
    
    library(brms)
    library(dplyr)

    thisEnv <- environment()
    fit <- NULL
    k <- NULL
    group_par <- gsub('[\\._]', '', objects)
    pars <- gsub('[\\._]', '', vars)

    if (type %in% 2:3) {
        if (type == 2)
            bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + (1 | ', objects, ') + (1 | ', judges, ')')))
        if (type == 3)
            bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + ', judges, ' + (1 | ', objects, ')')))
    } else {
        bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + (1 | ', objects, ')')))
    }
    
    if (length(vars) > 1)
        bform <- bform + set_rescor(rescor)

    me <- list(
        thisEnv = thisEnv,
        
        fit = function(df, ..., prior = NULL, rescor = F, cores = 4, k = NULL, k_fun = min, k_fun_args = list(na.rm=T)) {
            # Fits the multivariate Bayesian model using brms
            #
            # Args:
            #  df          Data frame in which to find the variables
            #  k_fun       Function to calculate the ICC*k values from the table of observations per
            #              variable and per object (the default, min, finds the minimum number of observations
            #              of each object). Can be overriden by setting k when calling icc()
            #  k_fun_args  List of named arguments passed on to k_fun.
            #  ...         Arguments passed to brms
            fit <- get('fit', thisEnv)
            if (! is.null(fit))
                message("Recompiling before refitting; use refit() to avoid recompilation")
            if (is.null(prior)) {
                m <- get('me', thisEnv)
                prior <- m$get_prior(df)
            }
            bform <- get('bform', thisEnv)
            vars <- get('vars', thisEnv)
            objects <- get('objects', thisEnv)
            result <- brm(bform, data = df, cores = cores, prior = prior, ...)
            assign('fit', result, envir = thisEnv)
            k <- do.call(apply, c(list(X = df %>% select_(.dots = c(objects, vars)) %>% group_by_(objects) %>% summarise_all(funs(sum(!is.na(.)))) %>% select(-1), MARGIN=2, FUN=k_fun), k_fun_args))
            assign('k', k, envir = thisEnv)
            return(result)
        },
        
        refit = function(..., cores = 4, overwrite = T) {
            # Refit the same model to the same data, e.g., with new parameters
            #
            # Args:
            #  overwrite  logical, whether or not to overwrite the internal fit object (default T)
            #
            # Returns:
            #  the new fit
            fit <- get('fit', thisEnv)
            result <- update(fit, cores = cores, ...)
            if (overwrite)
                assign('fit', result, envir = thisEnv)
            return(result)
        },
        
        get_fit = function() {
            # Get the existing fit
            fit <- get('fit', thisEnv)
            if (is.null(fit))
                stop("Run method 'fit' first")
            return(fit)
        },

        get_k = function() {
            # Get the divider
            k <- get('k', thisEnv)
            if (is.null(k))
                stop("Run method 'fit' first")
            return(k)
        },
        
        get_prior = function(df, sd_prior = NULL, priors=c('cauchy', 'default', 'normal', 'uniform', 'inv_gamma')) {
            # Set priors for the sd variables to sd_prior
            # or auto-choose a prior
            # Note that this leaves the priors for sigma alone, which are student t by default scaled around the approx. mean of the variable
            bform <- get('bform', thisEnv)
            priors <- match.arg(priors)
            p <- brms::get_prior(bform, df)
            if (! is.null(sd_prior)) {
                p[p$class == 'sd', 'prior'] <- sd_prior
            } else {
                if (priors == 'cauchy') {
                    p[p$class == 'sd', 'prior'] <- 'cauchy(0, 5)'
                } else if (priors == 'inv_gamma') {
                    p[p$class == 'sd', 'prior'] <- 'inv_gamma(1, 1)'
                } else if (priors == 'normal') {
                    p[p$class == 'sd', 'prior'] <- 'normal(0, 5)'
                } else if (priors == 'uniform') {
                    # Erase all default priors because brms then defaults to uniform on [0, ∞]
                    # which is recommended according to doi:10.1186/1471-2288-14-121
                    # however Gelman (2006) and others find that this over-estimates the variance components
                    p[p$class == 'sd', 'prior'] <- ''
                }
            }
            return(p)
        },

        get_bform = function() {
            return(get('bform', thisEnv))
        },

        set_bform = function(bform) {
            message('Use with extreme caution, YMMV')
            assign('bform', bform, envir = thisEnv)
        },

        # The functions below generate the calculations to estimate the ICCs
        #
        # Args:
        #  k   Optional vector of k values; default is to determine from the data
        #
        # Returns:
        #   a character vector of calculations
        
        icc1_formulae = function(k = NULL) {
            if (get('type', thisEnv) != 1)
                warning("Model is not set up for ICC1 but for ICC", get('type', thisEnv), "; use result with care")
            if (is.null(k))
                k <- get('k', thisEnv)
            vars <- get('vars', thisEnv)
            objects <- get('objects', thisEnv)
            pars <- get('pars', thisEnv)
            if (length(pars) > 1) {
                pars <- paste0('_', pars)
            } else {
                pars <- ''
            }
            formulae <- c(
                paste0('sd_', objects, '_', pars, '_Intercept^2 / (sd_', objects, '_', pars, '_Intercept^2 + sigma', pars, '^2)'),
                paste0('sd_', objects, '_', pars, '_Intercept^2 / (sd_', objects, '_', pars, '_Intercept^2 + sigma', pars, '^2 / ', k, ')')
            )
            names(formulae) <- paste(rep(vars, 2), rep(c('ICC1', 'ICC1k'), each=length(vars)))
            return(formulae)
        },
        
        icc2_formulae = function(k = NULL) {
            if (get('type', thisEnv) == 1)
                stop("Model is not set up for ICC2 but for ICC1; cannot generate hypotheses.")
            if (get('type', thisEnv) == 3)
                warning("Model is not set up for ICC2 but for ICC3; generated hypotheses cannot be tested in this model.")
            if (is.null(k))
                k <- get('k', thisEnv)
            vars <- get('vars', thisEnv)
            objects <- get('objects', thisEnv)
            judges <- get('judges', thisEnv)
            pars <- get('pars', thisEnv)
            if (length(pars) > 1) {
                pars <- paste0('_', pars)
            } else {
                pars <- ''
            }
            formulae <- c(
                paste0('sd_', objects, '_', pars, '_Intercept^2 / (sd_', objects, '_', pars, '_Intercept^2 + sd_', judges, '_', pars, '_Intercept^2 + sigma', pars, '^2)'),
                paste0('sd_', objects, '_', pars, '_Intercept^2 / (sd_', objects, '_', pars, '_Intercept^2 + sd_', judges, '_', pars, '_Intercept^2 + sigma', pars, '^2 / ', k, ')')
            )
            names(formulae) <- paste(rep(vars, 2), rep(c('ICC2', 'ICC2k'), each=length(vars)))
            return(formulae)
        },
        
        icc3_formulae = function(k = NULL) {
            # the calculation of ICC3 is identical to the correlation of ICC1
            # but the result is different because the model is different
            m <- get('me', thisEnv)
            formulae <- suppressWarnings(m$icc1_formulae(k = k))
            names(formulae) <- sub('ICC1', 'ICC3', names(formulae))
            return(formulae)
        },

        icc = function(test = "= 0", k = NULL) {
            # The functions below estimate the actual ICCs
            #
            # Args:
            #  test_values   The value(s) against which the calculated ICCs are compared
            #
            # Returns:
            #   a brmshypothesis object, that can be used to plot or to extract confidence intervals
            m <- get('me', thisEnv)
            fit <- m$get_fit()
            if (get('type', thisEnv) == 1) {
                fs <- m$icc1_formulae(k = k)
            } else if (get('type', thisEnv) == 2) {
                fs <- m$icc2_formulae(k = k)
            } else if (get('type', thisEnv) == 3) {
                fs <- m$icc3_formulae(k = k)
            }
            hyps <- paste(fs, test)
            names(hyps) <- names(fs)
            hypothesis(fit, hyps, class=NULL)
        },

        icc.tidyMCMC = function(..., test = "= 0", k = NULL) {
            # Return a tidy version of the iccs
            library(broom)
            m <- get('me', thisEnv)
            icc <- m$icc()
            bind_cols(data_frame(par = rownames(icc$hypothesis)), tidyMCMC(icc$samples, ...)) %>%
                select(-term)
        },

        getEnv = function() {
            return(get('thisEnv', thisEnv))
        },
        
        plot = function(conf.method = 'HPDinterval', conf.level = .95, regex_pars = '.*') {
            library(ggplot2)
            library(broom)
            library(dplyr)
            m <- get('me', thisEnv)
            icc <- m$icc()
            hyps <- rownames(icc$hypothesis)
            names(hyps) <- paste0('H', 1:length(hyps))
            tidyMCMC(icc$samples, conf.int = T, conf.method = conf.method, conf.level = conf.level) %>%
                filter(grepl(regex_pars, hyps[term])) %>%
                mutate(par = sub('\\sICC[123]k?$', '', hyps[term]), icc = sub('.*\\s(ICC[123]k?)$', '\\1', hyps[term])) %>%
                ggplot(aes(x = par, y = estimate)) +
                geom_violin(data = as_data_frame(icc$samples) %>%
                                gather(term, estimate) %>%
                                filter(grepl(regex_pars, hyps[term])) %>%
                                mutate(par = sub('\\sICC[123]k?$', '', hyps[term]), icc = sub('.*\\s(ICC[123]k?)$', '\\1', hyps[term])),
                            fill = 'white', color = 'lightblue', width=1) +
                geom_pointrange(aes(ymin = conf.low, ymax = conf.high), color='darkblue') +
                coord_flip() +
                facet_grid(~ icc) +
                labs(y = "ICC Posterior densities, medians,\nand 95% credible intervals", x = NULL) +
                scale_y_continuous(expand = c(0, .01), limits=c(0, 1), breaks = seq(0, 1, len=6)) +
                theme(panel.spacing = unit(1.5, 'lines'))
        }
    )

    assign('this', me, envir = thisEnv)
    class(me) <- append(class(me), 'bayesian_ICC')

    return(me)
}
