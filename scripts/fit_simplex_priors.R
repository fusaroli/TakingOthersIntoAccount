## function to fit a circular inference model through brms with the following modifications:
# - simplex reparametrisation of participant-level sigma
# - priors for the intercepts on the correct scale

## ## copy of brms to allow for manual tweaking of code
## brm_fit_code <- function (formula, data, code, family = gaussian(), prior = NULL, autocor = NULL, 
##                           cov_ranef = NULL, sample_prior = c("no", "yes", "only"), 
##                           sparse = FALSE, knots = NULL, stanvars = NULL, stan_funs = NULL, 
##                           fit = NA, save_ranef = TRUE, save_mevars = FALSE, save_all_pars = FALSE, 
##                           inits = "random", chains = 4, iter = 2000, warmup = floor(iter/2), 
##                           thin = 1, cores = getOption("mc.cores", 1L), control = NULL, 
##                           algorithm = c("sampling", "meanfield", "fullrank"), future = getOption("future", 
##                                                                                                  FALSE), silent = TRUE, seed = NA, save_model = NULL, 
##                           save_dso = TRUE, ...) 
## {
##     ## from some old brms version (~2018)
##     dots <- list(...)
##     autocor <- brms:::check_autocor(autocor)
##     algorithm <- match.arg(algorithm)
##     testmode <- dots$testmode
##     dots$testmode <- NULL
##     if (is(fit, "brmsfit")) {
##         x <- fit
##         sdata <- standata(x)
##         x$fit <- rstan::get_stanmodel(x$fit)
##     }
##     else {
##         formula <- brms:::validate_formula(formula, data = data, family = family, 
##                                     autocor = autocor)
##         family <- brms:::get_element(formula, "family")
##         autocor <- brms:::get_element(formula, "autocor")
##         bterms <- parse_bf(formula)
##         if (is.null(dots$data.name)) {
##             data.name <- substr(collapse(deparse(substitute(data))), 
##                                 1, 50)
##         }
##         else {
##             data.name <- dots$data.name
##             dots$data.name <- NULL
##         }
##         data <- brms:::update_data(data, bterms = bterms)
##         prior <- brms:::check_prior(prior, formula, data = data, sample_prior = sample_prior, 
##                              warn = TRUE)
##         x <- brms:::brmsfit(formula = formula, family = family, data = data, 
##                      data.name = data.name, prior = prior, autocor = autocor, 
##                      cov_ranef = cov_ranef, stanvars = stanvars, stan_funs = stan_funs, 
##                      algorithm = algorithm)
##         x$ranef <- brms:::tidy_ranef(bterms, data = x$data)
##         x$exclude <- brms:::exclude_pars(bterms, data = x$data, ranef = x$ranef, 
##                                   save_ranef = save_ranef, save_mevars = save_mevars, 
##                                   save_all_pars = save_all_pars)
##         # x$model <- make_stancode(formula, data = data, prior = prior, 
##         #                          sparse = sparse, cov_ranef = cov_ranef, sample_prior = sample_prior, 
##         #                          knots = knots, stanvars = stanvars, stan_funs = stan_funs, 
##         #                          save_model = save_model, brm_call = TRUE)
##         sdata <- make_standata(formula, data = data, prior = prior, 
##                                cov_ranef = cov_ranef, sample_prior = sample_prior, 
##                                knots = knots, stanvars = stanvars)
##         message("Compiling the C++ model")
##         x$fit <- brms:::eval_silent(rstan::stan_model(model_code = code, 
##                                                       save_dso = save_dso))
##         x$model <- code
##     }
##     if (is.character(inits) && !inits %in% c("random", "0")) {
##         inits <- get(inits, mode = "function", envir = parent.frame())
##     }
##     args <- nlist(object = x$fit, data = sdata, pars = x$exclude, 
##                   include = FALSE, algorithm, iter, seed)
##     args[names(dots)] <- dots
##     message("Start sampling")
##     if (args$algorithm == "sampling") {
##         args$algorithm <- NULL
##         args <- c(args, nlist(init = inits, warmup, thin, control, 
##                               show_messages = !silent))
##         if (future) {
##             brms:::require_package("future")
##             if (cores > 1L) {
##                 brms:::warning2("Argument 'cores' is ignored when using 'future'.")
##             }
##             args$chains <- 1L
##             futures <- fits <- vector("list", chains)
##             for (i in seq_len(chains)) {
##                 args$chain_id <- i
##                 if (is.list(inits)) {
##                     args$init <- inits[i]
##                 }
##                 futures[[i]] <- future::future(do.call(rstan::sampling, 
##                                                        args), packages = "rstan")
##             }
##             for (i in seq_len(chains)) {
##                 fits[[i]] <- future::value(futures[[i]])
##             }
##             x$fit <- rstan::sflist2stanfit(fits)
##             rm(futures, fits)
##         }
##         else {
##             args <- c(args, nlist(chains, cores))
##             x$fit <- do.call(rstan::sampling, args)
##         }
##     }
##     else {
##         x$fit <- do.call(rstan::vb, args)
##     }
##     if (!isTRUE(testmode)) {
##         x <- brms:::rename_pars(x)
##     }
##     x
## }

fit_brm = function(...) {
    # extracts the stancode from the model,
    # then replaces the definition of the standard deviation of the random effects with
    # one common scale and a simplex to rule them
    
    
    code = make_stancode(...)
    
    which_random = 1 # TODO more than one random effects group
    
    #################
    # simplex parametrization
    for (i in which_random) {
        code = code %>%
            

        ##
        # replace the parameter sd_1 with sd_1_simplex and sd_1_scale
            str_replace("vector<lower=0>\\[(M_\\d)\\] sd_(\\d);  // group-level standard deviations",
                        "simplex[\\1] sd_\\2_simplex; // distribution of group-level standard deviations;
                         real<lower=0> sd_\\2_scale; // scale of group-level standard deviations") %>%

        ##
        #add sd_1 as a transformed parameter
            str_replace("transformed parameters \\{",
                        glue("transformed parameters {
                              vector<lower=0>[M_@i!] sd_@i! = sd_@i!_simplex * sd_@i!_scale;", .open = "@", .close = "!")) %>%

        ##
        # replace the prior for sd_1 with a prior for sd_1_scale
            str_replace_all(glue("\\(sd_@i! \\|",   .open="@", .close="!"),
                        glue("(sd_@i!_scale |", .open="@", .close="!"))


        

    #######################
    # priors on correct scale
        
    }
    w_params = str_match_all(code, "vector\\[.*?\\] b_(w\\w+)")[[1]][,2]
    a_params = str_match_all(code, "vector\\[.*?\\] b_(a\\w+)")[[1]][,2]
    for (w in w_params) {
        n_coef = length(str_match_all(code, glue("b_{w}\\[[23456789]\\d*\\]"))[[1]])
        code <- code %>%
            # declare the intercept and the coefficient vector
            str_replace("parameters \\{",
                        glue("parameters {{\nreal<lower=.5, upper=1> b_{w}_Intercept_unscaled;")) %>%
            
            # change name of which params to put priors on
            str_replace(glue("{w}(\\[1\\])? \\|"),
                        glue("{w}_Intercept_unscaled |"))  %>%
            
            str_replace(glue("vector\\[K_{w}\\] b_{w};  // population-level effects"),
                        "") %>%
            
            str_replace("transformed parameters \\{",
                        glue("transformed parameters {{\n  vector[K_{w}] b_{w};")) %>%
            
            str_replace("\\}.*?\nmodel \\{",
                        glue("  b_{w}[1] = logit((b_{w}_Intercept_unscaled-.5 )*2);\n}\nmodel {{"))
        
        if (n_coef > 0) {
            # more than just intercept
            code = code %>%
                str_replace("parameters \\{",
                            glue("parameters {{\nvector[K_{w}-1] b_{w}_coef;")) %>%
                
                str_replace(glue("(target \\+= \\w+\\()b_{w}\\[[23456789]\\d*\\] \\| (\\w+), (\\d+(\\.\\d+)?)\\);"),
                            glue("\\1b_{w}_coef | \\2, \\3);"))  %>%
                
                str_replace_all(glue("(target \\+= \\w+\\()b_{w}\\[[23456789]\\d*\\] \\| (\\w+), (\\d+(\\.\\d+)?)\\);"),
                                glue(""))  %>%
                
                str_replace("\\}.*?\nmodel \\{",
                            glue("b_{w}[2:K_{w}] = b_{w}_coef;\n}}\nmodel {{"))
            
        }
    }
    
    for (a in a_params) {
        n_coef = length(str_match_all(code, glue("b_{a}\\[[23456789]\\d*\\]"))[[1]])
        code <- code %>%
            # declare the intercept and the coefficient vector
            str_replace("parameters \\{",
                        glue("parameters {{\nreal<lower=0> b_{a}_Intercept_unscaled;")) %>%
            
            # change name of which params to put priors on
            str_replace(glue("{a}(\\[1\\])? \\|"),
                        glue("{a}_Intercept_unscaled |"))  %>%
            
            str_replace(glue("vector\\[K_{a}\\] b_{a};  // population-level effects"),
                        "") %>%
            
            str_replace("transformed parameters \\{",
                        glue("transformed parameters {{\nvector[K_{a}] b_{a};")) %>%
            
            str_replace("\\}.*?\nmodel \\{",
                        glue("b_{a}[1] = log(b_{a}_Intercept_unscaled);\n}\nmodel {{"))
        
        if (n_coef > 0) {
            # more than just intercept
            code = code %>%
                str_replace("parameters \\{",
                            glue("parameters {{\nvector[K_{a}-1] b_{a}_coef;")) %>%
                
                str_replace(glue("(target \\+= \\w+\\()b_{a}\\[[23456789]\\d*\\] \\| (\\w+), (\\d+(\\.\\d+)?)\\);"),
                            glue("\\1b_{a}_coef | \\2, \\3);"))  %>%
                
                str_replace_all(glue("(target \\+= \\w+\\()b_{a}\\[[23456789]\\d*\\] \\| (\\w+), (\\d+(\\.\\d+)?)\\);"),
                                glue(""))  %>%
                
                str_replace("\\}.*?\nmodel \\{",
                            glue("b_{a}[2:K_{a}] = b_{a}_coef;\n}}\nmodel {{"))
            
        }
    }

    
    
    
    cat(code, file = "models/tmp_brms_model.stan")
    
    dots = list(...)
    dots$code = code
    do.call(brm_fit_from_code, dots)
}


brm_fit_from_code <- function (formula, data, code, family = gaussian(), prior = NULL, autocor = NULL, 
    data2 = NULL, cov_ranef = NULL, sample_prior = c("no", "yes", 
        "only"), sparse = NULL, knots = NULL, stanvars = NULL, 
    stan_funs = NULL, fit = NA, save_ranef = TRUE, save_mevars = FALSE, 
    save_all_pars = FALSE, inits = "random", chains = 2, iter = 2000, 
    warmup = floor(iter/2), thin = 1, cores = 2, control = list(
        adapt_delta = .999,
        max_treedepth = 20), algorithm = c("sampling", "meanfield", 
        "fullrank"), future = getOption("future", FALSE), silent = TRUE, 
    seed = NA, save_model = NULL, stan_model_args = list(), save_dso = TRUE, 
    file = NULL, ...) 
{
    ## from brms version 2.12.0
    if (!is.null(file)) {
        x <- read_brmsfit(file)
        if (!is.null(x)) {
            return(x)
        }
    }
    dots <- list(...)
    testmode <- isTRUE(dots$testmode)
    dots$testmode <- NULL
    algorithm <- match.arg(algorithm)
    silent <- brms:::as_one_logical(silent)
    iter <- brms:::as_one_numeric(iter)
    warmup <- brms:::as_one_numeric(warmup)
    thin <- brms:::as_one_numeric(thin)
    chains <- brms:::as_one_numeric(chains)
    cores <- brms:::as_one_numeric(cores)
    future <- brms:::as_one_logical(future) && chains > 0L
    seed <- brms:::as_one_numeric(seed, allow_na = TRUE)
    if (is.character(inits) && !inits %in% c("random", "0")) {
        inits <- get(inits, mode = "function", envir = parent.frame())
    }
    if (is.brmsfit(fit)) {
        x <- fit
        icnames <- c("loo", "waic", "kfold", "R2", "marglik")
        x[icnames] <- list(NULL)
        sdata <- standata(x)
        x$fit <- rstan::get_stanmodel(x$fit)
    }
    else {
        formula <- brms:::validate_formula(formula, data = data, family = family, 
            autocor = autocor, sparse = sparse)
        family <- brms:::get_element(formula, "family")
        bterms <- parse_bf(formula)
        data.name <- brms:::substitute_name(data)
        data <- brms:::validate_data(data, bterms = bterms)
        data2 <- brms:::validate_data2(data2, bterms = bterms, brms:::get_data2_autocor(formula))
        prior <- brms:::check_prior(prior, formula = formula, data = data, 
            sample_prior = sample_prior, warn = FALSE)
        x <- brms:::brmsfit(formula = formula, family = family, data = data, 
            data.name = data.name, data2 = data2, prior = prior, 
            cov_ranef = cov_ranef, stanvars = stanvars, stan_funs = stan_funs, 
            algorithm = algorithm)
        x$ranef <- brms:::tidy_ranef(bterms, data = x$data)
        x$exclude <- brms:::exclude_pars(x, save_ranef = save_ranef, 
            save_mevars = save_mevars, save_all_pars = save_all_pars)
        ## x$model <- make_stancode(formula, data = data, prior = prior, 
        ##     cov_ranef = cov_ranef, sample_prior = sample_prior, 
        ##     knots = knots, stanvars = stanvars, stan_funs = stan_funs, 
        ##     save_model = save_model) ## malte edit
        x$model <- code ## malte edit
        sdata <- brms:::make_standata(formula, data = data, prior = prior, 
            data2 = data2, cov_ranef = cov_ranef, sample_prior = sample_prior, 
            knots = knots, stanvars = stanvars)
        stopifnot(is.list(stan_model_args))
        silence_stan_model <- !length(stan_model_args)
        stan_model_args$model_code <- x$model
        if (!isTRUE(save_dso)) {
            warning2("'save_dso' is deprecated. Please use 'stan_model_args'.")
            stan_model_args$save_dso <- save_dso
        }
        message("Compiling the C++ model")
        x$fit <- brms:::eval_silent(do_call(rstan::stan_model, stan_model_args), 
            silent = silence_stan_model, type = "message")
    }
    args <- nlist(object = x$fit, data = sdata, pars = x$exclude, 
        include = FALSE, algorithm, iter, seed)
    args[names(dots)] <- dots
    message("Start sampling")
    if (args$algorithm == "sampling") {
        args$algorithm <- NULL
        brms:::`c<-`(args, nlist(init = inits, warmup, thin, control, 
            show_messages = !silent))
        if (future) {
            if (cores > 1L) {
                brms:::warning2("Argument 'cores' is ignored when using 'future'.")
            }
            args$chains <- 1L
            futures <- fits <- vector("list", chains)
            for (i in seq_len(chains)) {
                args$chain_id <- i
                if (is.list(inits)) {
                  args$init <- inits[i]
                }
                futures[[i]] <- future::future(brms::do_call(rstan::sampling, 
                  args), packages = "rstan")
            }
            for (i in seq_len(chains)) {
                fits[[i]] <- future::value(futures[[i]])
            }
            x$fit <- rstan::sflist2stanfit(fits)
            rm(futures, fits)
        }
        else {
            brms:::`c<-`(args, nlist(chains, cores))
            x$fit <- do_call(rstan::sampling, args)
        }
    }
    else {
        x$fit <- do_call(rstan::vb, args)
    }
    if (!testmode) {
        x <- brms:::rename_pars(x)
    }
    if (!is.null(file)) {
        brms:::write_brmsfit(x, file)
    }
    x
}
