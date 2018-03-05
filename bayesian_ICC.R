bayesian_ICC <- function(vars, objects, judges = NA, which = 1, rescor = F) {
    # Class to fit six ICCs (ICC1, ICC1k, ICC2, ICC2k, ICC3, ICC3k) on multiple variables simultaneously in a Bayesian framework with ease
    #
    # Args:
    #  vars      vector of characters, specifying the variable names of the variables to aggregate
    #  objects   character, specifying the objects having been rated (e.g., the number of groups or video clips)
    #  judges    character, optional specification for specification of judges (required to calculate ICC2 and ICC3)
    #  rescor    Whether or not to model response correlation if > 1 var
    #
    # Returns:
    #  A class with multiple methods
    
    if (! is.character(vars) | ! is.character(objects) | length(objects) > 1)
        stop("Vars should be a character vector; objects should be a character variable")
    if (! is.na(judges) & (!is.character(judges) | length(judges) > 1))
        stop("Judges should be NA or a character variable")
    if (length(which) != 1 | ! which %in% 1:3)
        stop("Which should be 1, 2, or 3 for ICC1, ICC2, and ICC3, respectively")
    if (which %in% 2:3 & is.na(judges))
        stop("Specify judges for ICC2 or ICC3")
    if (which == 1 & !is.na(judges))
        warning("Judges argument will be ignored when calculating ICC1")
    
    library(brms)
    library(dplyr)

    thisEnv <- environment()
    fit <- NULL
    dividers <- NULL
    group_par <- gsub('[\\._]', '', objects)
    pars <- gsub('[\\._]', '', vars)

    if (which %in% 2:3) {
        if (which == 2)
            bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + (1 | ', objects, ') + (1 | ', judges, ')')))
        if (which == 3)
            bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + judges + (1 | ', objects, ')')))
    } else {
        bform <- bf(as.formula(paste0('cbind(', paste(vars, collapse=', '), ') ~ 1 + (1 | ', objects, ')')))
    }
    
    if (length(vars) > 1)
        bform <- bform + set_rescor(rescor)

    me <- list(
        thisEnv = thisEnv,
        
        fit = function(df, ..., rescor = F, cores = 4, k = NA, k_fun = min, k_fun_args = list(na.rm=T)) {
            # Fits the multivariate Bayesian model using brms
            #
            # Args:
            #  df          Data frame in which to find the variables
            #  k           The k value used to calculate the ICC*k values. This is normally the number of observations of each object (e.g., the number of group members). If NA, k_fun and k_fun_args will be used to find this value.
            #  k_fun       Function to run on the table of observations per variable and per object (the default, min, finds the minimum number of observations of each object)
            #  k_fun_args  List of named arguments passed on to k_fun.
            #  ...         Arguments passed to brms
            fit <- get('fit', thisEnv)
            if (! is.null(fit))
                message("Refitting; previous fit will be overwritten")
            bform <- get('bform', thisEnv)
            vars <- get('vars', thisEnv)
            objects <- get('objects', thisEnv)
            result <- brm(bform, data = df, cores = cores, ...)
            assign('fit', result, envir = thisEnv)
            dividers <- do.call(apply, c(list(X = df %>% select_(.dots = c(objects, vars)) %>% group_by_(objects) %>% summarise_all(funs(sum(!is.na(.)))) %>% select(-1), MARGIN=2, FUN=k_fun), k_fun_args))
            assign('dividers', dividers, envir = thisEnv)
            return(result)
        },
        
        refit = function(overwrite = T, ..., cores = 4) {
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

        # The functions below generate the hypotheses to test ICC1s
        #
        # Args:
        #  test_value   The value against which the calculated ICC is compared
        #
        # Returns:
        #   a character vector of hypotheses
        
        icc1_hypotheses = function(test_value = .01) {
            m <- get('me', thisEnv)
            m$icc1k_hypotheses(test_value = test_value, dividers = 1)
        },
        
        icc1k_hypotheses = function(test_value = .01, dividers = NA) {
            if (is.na(dividers))
                dividers <- get('dividers', thisEnv)
            vars <- get('vars', thisEnv)
            objects <- get('objects', thisEnv)
            if (length(pars) > 1) {
                pars <- paste0('_', get('pars', thisEnv))
            } else {
                pars <- ''
            }
            return(hyps.icc1k <- paste0('sd_', objects, '_', pars, '_Intercept^2 / (sd_', objects, '_', pars, '_Intercept^2 + sigma', pars, '^2 / ', dividers, ') > ', test_value))
        },
        
        icc2_hypothesis = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc2k_hypothesis = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc3_hypothesis = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc3k_hypothesis = function(test_value = .01) {
            stop("Not implemented yet")
        },
        
        # The functions below estimate the actual ICCs
        #
        # Args:
        #  test_value   The value against which the calculated ICC is compared
        #
        # Returns:
        #   a brmshypothesis object, that can be used to plot or to extract confidence intervals
        icc1 = function(test_value = .01) {
            m <- get('me', thisEnv)
            fit <- m$get_fit()
            hyp <- m$icc1_hypotheses(test_value = test_value)
            return(hypothesis(fit, hyp, class = NULL))
        },
        
        icc1k = function(test_value = .01) {
            m <- get('me', thisEnv)
            fit <- m$get_fit()
            hyp <- m$icc1k_hypotheses(test_value = test_value)
            return(hypothesis(fit, hyp, class = NULL))
        },
        
        icc2 = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc2k = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc3 = function(test_value = .01) {
            stop("Not implemented yet")
        },
        icc3k = function(test_value = .01) {
            stop("Not implemented yet")
        },
        
        getEnv = function() {
            return(get('thisEnv', thisEnv))
        }
    )
    assign('this', me, envir = thisEnv)
    class(me) <- append(class(me), 'bayesian_ICC')
    return(me)
}
