library(tidyverse)
library(brms)
# library(tidybayes)
library(bayesplot)
library(boot)
library(rlang)
### Usage: 
# plot_ppc(circular_m, "prior")
# plot_ppc(circular_m, "posterior")
# 
# plot_ppc(circular_d_m, "prior", diagnosis = TRUE)
# plot_ppc(circular_d_m, "posterior", diagnosis = TRUE)

table_of_social_evidence = full_data %>%
    distinct(lOthers1, lOthers2, lOthers3, lOthers4, lOthers,
             lOthers1NoC, Confidence1, lOthers2NoC, Confidence2, lOthers3NoC, Confidence3, lOthers4NoC, Confidence4,
             lOthersM) %>%
    rowwise() %>%
    mutate(conf_mean = mean(c(Confidence1, Confidence2, Confidence3, Confidence4))) %>%
    ungroup()



combine_others <- function(dat) {
    dat %>%
        left_join(table_of_social_evidence)
}

add_prior_samples  <- function(newdata, model, diagnosis, re_prior, re_name,...) {
    # re_name = as_string(re_name)
    s = posterior_samples(model, pars = "prior_b")
    n_params = length(names(s))
    print(names(s))
    names(s) = str_sub(names(s), 9)
    
    # manually sample from the random effects prior
    # student_t(0, 3)
    # lkj(8)
    if (is.na(re_prior[[1]])) re_prior = c(3, 8)
    # sd = rt(n_params, pluck(re_prior, 1)) / n_params
    n_participants = length(unique(select(newdata, !!re_name) %>% pluck(1)))
    sd = rnorm(n_params, 0, pluck(re_prior, 1)) / n_params
    omega = rethinking::rlkjcorr(n_participants, n_params, eta = pluck(re_prior, 2)) 
    Sigma = purrr::map(1:n_participants, ~ diag(sd) %*% omega[.,,] %*% diag(sd))
    # Sigma = diag(sd) %*% omega %*% diag(sd)
    re = map(1:length(unique(newdata[as_string(re_name)]) %>% pluck(1)),
             # for each participant ...
             function(p) {
                 result = MASS::mvrnorm(nsamples(model),
                                        rep(0, n_params),
                                        Sigma[[p]]) %>%
                     as.data.frame()
                 names(result) = str_c(names(s), "_re")
                 return(result)
             })
    
    # print(str(re))
    # print(re)
    # names(re) = unique(newdata[as_string(re_name)])
    s = mutate(s, .iteration = 1:n())
    # 
    
    
    dat = mutate(newdata, s = list(s), re = list(re),
                 .row = 1:n()) %>% 
        mutate(re = map2(re, !!re_name, pluck)) %>%
        unnest()
    
    if (diagnosis) {
        dat = dat %>% 
            group_by(.row) %>%
            mutate(
                bias_Intercept_re = bias_re,
                bias_Diagnosis_re = sample(bias_re),
                wSelf_Intercept_re = logit_scaled(wSelf_re, .5, 1),
                wSelf_Diagnosis_re = sample(wSelf_re),
                wOthers_Intercept_re = logit_scaled(wOthers_re, .5, 1),
                wOthers_Diagnosis_re = sample(wOthers_re),
                
                bias_Intercept = bias,
                bias_Diagnosis = sample(bias),
                wSelf_Intercept = logit_scaled(wSelf, .5, 1),
                wSelf_Diagnosis = sample(wSelf),
                wOthers_Intercept = logit_scaled(wOthers, .5, 1),
                wOthers_Diagnosis = sample(wOthers))
        if ("aOthers" %in% names(s)) {
            dat = mutate(dat,
                 aOthers_Intercept_re = log(aOthers_re),
                 aOthers_Diagnosis_re = sample(aOthers_re),
                 aSelf_Intercept_re = log(aSelf_re),
                 aSelf_Diagnosis_re = sample(aSelf_re),
                 
                 aOthers_Intercept = log(aOthers),
                 aOthers_Diagnosis = sample(aOthers),
                 aSelf_Intercept = log(aSelf),
                 aSelf_Diagnosis = sample(aSelf))
        }
                         
    } else {
        dat = dat %>% 
                group_by(.row) %>%
                mutate(
                    bias_Intercept_re = bias_re,
                    wSelf_Intercept_re = wSelf_re,
                    wOthers_Intercept_re = wOthers_re,
                    
                    bias_Intercept = bias,
                    wSelf_Intercept = wSelf,
                    wOthers_Intercept = wOthers)
            if ("aOthers" %in% names(s)) {
                dat = mutate(dat,
                             aOthers_Intercept_re = aOthers_re,
                             aSelf_Intercept_re = aSelf_re,
                             
                             aOthers_Intercept = aOthers,
                             aSelf_Intercept = aSelf)
                
            }
        }
    
    ungroup(dat)
}   

add_posterior_samples <- function(newdata, model, diagnosis, p_names, re_name, ...) {
    s = posterior_samples(model, pars = "^b_")
    names(s) = str_sub(names(s), 3)
    s = mutate(s, .iteration = 1:n())
    
    # extract random effects
    # dims: iteration, participant, param
    re = ranef(model, summary = FALSE)[[1]] %>%
        apply(2, list) %>%
        map(function(x) as.data.frame(x) %>%
            rename_all(function(name) {str_c(name, "_re")}))
    # 
    # print(s)
    # print(re)
    
    # now a list of participant data frames, each with iter rows and nparam cols
    dat = mutate(newdata, s = list(s), re = list(re)) %>%
        # print()
        # select the correct random effect for a given participant
        mutate(re = map2(re, !!re_name, function(r, p) {
            pluck(r, as.character(pluck(p_names, p)))})) %>% 
        unnest() 
    
    # View(dat)
    # print(str(dat))
    return(dat)
}

integrate_diagnosis <- function(dat, diagnosis) {
    if (diagnosis) {
        try({
            dat = dat %>%
                mutate(aSelf = aSelf_Intercept + aSelf_Diagnosis * Diagnosis,
                       aOthers = aOthers_Intercept + aOthers_Diagnosis * Diagnosis)
            
        }, silent = TRUE)
        
        try({
            
            dat = dat %>%
                mutate(bConfidence = bConfidence_Intercept + bConfidence_Diagnosis * Diagnosis)
            
        }, silent = TRUE)
        
        dat = dat %>%
            mutate(bias = bias_Intercept + bias_Diagnosis * Diagnosis,
                   wSelf = wSelf_Intercept + wSelf_Diagnosis * Diagnosis,
                   wOthers = wOthers_Intercept + wOthers_Diagnosis * Diagnosis)
        
    } else {
        try({
            dat = dat %>%
                mutate(aSelf = aSelf_Intercept,
                       aOthers = aOthers_Intercept)
            
        }, silent = TRUE)
        try({
            dat = dat %>%
                mutate(bConfidence = bConfidence_Intercept)
            
        }, silent = TRUE)
        dat = dat %>%
            mutate(bias = bias_Intercept,
                   wSelf = wSelf_Intercept ,
                   wOthers = wOthers_Intercept)
        
    }
    
    return(dat)
}

add_bias <- function(dat) {
    if (!"bias" %in% names(dat)) {
        dat$bias = 0
    }
    return(dat)
}

add_circular_model <- function(dat, confidence) {
    if (confidence) {
    try({
    dat = mutate(dat,
           I = F(aSelf, lSelf, wSelf) +
               F(aOthers, lOthers1, wOthers + bConfidence * Confidence1) +
               F(aOthers, lOthers2, wOthers + bConfidence * Confidence2) +
               F(aOthers, lOthers3, wOthers + bConfidence * Confidence3) +
               F(aOthers, lOthers4, wOthers + bConfidence * Confidence4 ),
           circular = bias + 
               F(0, lSelf + I, wSelf) +
               F(0, lOthers1 + I, wOthers + bConfidence * Confidence1) +
               F(0, lOthers2 + I, wOthers + bConfidence * Confidence2) +
               F(0, lOthers3 + I, wOthers + bConfidence * Confidence3) +
               F(0, lOthers4 + I, wOthers + bConfidence * Confidence4),
           circular = inv.logit(circular)) %>%
        select(-I)
    }, silent = TRUE)
    } else {
     # no confidence   
        try({
            dat = mutate(dat,
                         I = F(aSelf, lSelf, wSelf) +
                             F(aOthers, lOthers1, wOthers) +
                             F(aOthers, lOthers2, wOthers) +
                             F(aOthers, lOthers3, wOthers) +
                             F(aOthers, lOthers4, wOthers ),
                         circular = bias + 
                             F(0, lSelf + I, wSelf) +
                             F(0, lOthers1 + I, wOthers) +
                             F(0, lOthers2 + I, wOthers) +
                             F(0, lOthers3 + I, wOthers) +
                             F(0, lOthers4 + I, wOthers),
                         circular = inv.logit(circular)) %>%
                select(-I)
        }, silent = TRUE)
        
    }
    
    return(dat)
}

add_naive_bayes <- function(dat) {
    mutate(dat, naive = inv.logit(bias + lSelf + lOthers1 + lOthers2 + lOthers3 + lOthers4))
}

add_weighted_bayes <- function(dat, diagnosis, confidence) {
    if (confidence) {
        dat =  dat %>%
            mutate(weighted = bias + 
                       F(0, lSelf, wSelf) +
                       F(0, lOthers1NoC, wOthers + bConfidence * Confidence1) +
                       F(0, lOthers2NoC, wOthers + bConfidence * Confidence2) +
                       F(0, lOthers3NoC, wOthers + bConfidence * Confidence3) +
                       F(0, lOthers4NoC, wOthers + bConfidence * Confidence4),
                   weighted = inv.logit(weighted))
        
    } else {
        dat =  dat %>%
            mutate(weighted = bias + 
                       F(0, lSelf, wSelf) +
                       F(0, lOthers1, wOthers) +
                       F(0, lOthers2, wOthers) +
                       F(0, lOthers3, wOthers) +
                       F(0, lOthers4, wOthers),
                   weighted = inv.logit(weighted))
    }
    return(dat)
    
}

round2 = function(...) as.character(round(as.numeric(...), digits = 2))

plot_ppc <- function(model, chose = circular, which = "posterior", diagnosis = FALSE,
                     random_effects = TRUE, by_participant = FALSE, participants = NA, re_prior = NA,
                     re_name = "Participant", confidence = FALSE, colour_name = "lOthersM",
                     plot_naive = TRUE) {
    
    which = match.arg(which, c("posterior", "prior"))
    re_name = sym(re_name)
    # colour_name = sym(colour_name)
    
    if (which == "posterior") {
        sample_function = add_posterior_samples
    } else {
        sample_function = add_prior_samples
    }
    if (diagnosis) {
        Diagnosis = c(0,1)
    } else {
        Diagnosis = 0
    }
    
    # use participant no 1,2,3 rather than 101, 103, 207 etc
    if (is.na(participants[[1]])) participants = unique(model$data[re_name])
    p_names = participants
    participants = 1:length(participants)
    # print(participants)
    
    
    # start building the data grid
    pl = modelr::data_grid(model$data,
                           # lSelf, lOthers1, lOthers2 = 0, lOthers3 =0, lOthers4 = 0,
                           lSelf = modelr::seq_range(lSelf, 20), table_of_social_evidence,
                           Diagnosis, 
                           Participant = participants) %>%
        dplyr::rename(!!re_name := Participant) %>%
        sample_function(model, diagnosis = diagnosis, re_prior = re_prior, p_names = p_names, re_name = re_name)
    
    
    if (random_effects) {
        has_re = str_match_all(names(pl), "(.*)_re$") %>%
            map(~pluck(., 2))
        
        for (param in has_re) {
            if (!is_null(param)) {
                # for every parameter (column) that also has a corresponding _re
                param_re = str_c(param, "_re")
                pl <- pl %>%
                    # add the _re column
                    mutate(!! param := !! ensym(param) + !! ensym(param_re))
            }
        }
    }
    
    pl = pl %>%
        integrate_diagnosis(diagnosis) %>%
        add_bias() %>%
        add_circular_model(confidence)    
    
    # prepare actual data
    actual_data = model$data %>%
        filter(!!re_name %in% p_names) %>%
        combine_others() %>%
        mutate(!!sym(colour_name) := factor(!!sym(colour_name))) %>%
        group_by(lSelf, !!sym(colour_name)) %>%
        mutate(n = n())
    
        
    
    pl = pl  %>%
        add_naive_bayes() %>%
        add_weighted_bayes(diagnosis, confidence) %>%
        # combine_others() %>%
        mutate(!!sym(colour_name) := factor(!!sym(colour_name)),
               Diagnosis = factor(Diagnosis, levels = c(0,1), labels = c("Control", "Diagnosis"))) %>%
        
        # rename participants to original names
        mutate(!!re_name := map_int(!!re_name, ~pluck(p_names, .))) %>%
        # plot predictions
        ggplot(aes_string("lSelf", quo_name(enquo(chose)), colour = colour_name, group = colour_name, fill = colour_name)) +
        geom_line(aes(linetype = "Model predictions"), stat = "summary", fun.y = mean) +
        geom_ribbon(stat = "summary", fun.data = mean_qi, alpha = .2, colour = NA) +
                        
        # uncomment the next two lines to include the naive bayes model
        # geom_ribbon(aes(y = naive), stat = "summary", fun.data = mean_qi, alpha = .1, colour = NA) +
        
            
            
        labs(subtitle = str_c(which, " predictive check (", ifelse(random_effects, "with", "without"), " random effects)"),
             title = str_c("Model: ", deparse(substitute(model))),
             y = "probability of chosing red",
             x = "log odds of non-social evidence",
             colour = "log odds of social evidence",
             fill = "log odds of social evidence") +
        coord_cartesian(xlim = c(-2.2, 2.2)) +
        scale_y_continuous(labels = scales::percent) +
        guides(linetype = guide_legend(title = element_blank())) +
        scale_color_brewer(palette = "RdYlBu", labels = round2, direction = -1) +
        scale_fill_brewer(palette = "RdYlBu", labels = round2, direction = -1)
    
    if (which == "posterior") {
        pl = pl + 
        # plot data
            geom_jitter(aes(y = ChoiceRed, group = str_c(lSelf, !!sym(colour_name))),
                        stat = "summary", fun.y = mean, data = actual_data,
                        width = .12, height = .01)
            # geom_pointrange(aes(y = ChoiceRed), stat = "summary", fun.data = mean_se, data = actual_data)
    }
    
    if (plot_naive) {
        pl = pl + 
            geom_line(aes(y = naive, linetype = "Naive bayes model"), stat="summary", fun.y=mean, alpha = .6)
    }
    
        # pl = case_when(diagnosis & by_participant ~ pl + facet_grid(by_participant ~ Diagnosis),
        #                 diagnosis ~ pl+facet_grid(~ Diagnosis),
        #           by_participant ~ pl+facet_wrap(~ Participant),
        #           TRUE ~ pl)
        # 
        # print(pl)
    

    if (diagnosis & by_participant) {
        pl + facet_grid(vars(!!re_name), vars(Diagnosis))
    } else if (diagnosis) {
        pl + facet_grid(~Diagnosis)
    } else if (by_participant) {
        pl + facet_wrap(vars(!!re_name))
    } else {
        pl
    }
    
}




#load("models/Arndis/confidence/weighted_confidence_diagnosis_m")
#ps = sample(unique(weighted_confidence_diagnosis_m$data$SubjectPair), 4)

#plot_ppc(weighted_confidence_diagnosis_m,
#          weighted,
#          "posterior",
#          diagnosis = TRUE,
#          random_effects = TRUE,
#          by_participant = TRUE,
#          participants = ps)
#ggsave("ppc_posterior_weighted_test1.png", width = 10, height = 9)
# 
# plot_ppc(weighted_m2,
#          weighted,
#          "posterior",
#          diagnosis = FALSE,
#          random_effects = TRUE,
#          by_participant = TRUE,
#          participants = ps)
# ggsave("ppc_posterior_weighted_test2.png", width = 10, height = 9)


# load("models/Arndis/circular_confidence_diagnosis_m")
# ps = sample(unique(circular_m$data$SubjectPair), 4)
# plot_ppc(circular_m,
#          circular,
#          "prior",
#          random_effects = TRUE,
#          by_participant = TRUE,
#          participants = ps)
# ggsave("ppc_posterior_circular_test.png", width = 10, height = 9)

# plot_ppc(circular_m,
#          circular,
#          "posterior",
#          random_effects = TRUE,
#          by_participant = TRUE,
#          participants = ps)
# ggsave("ppc_posterior_circular_test_re.png", width = 10, height = 9)
# 
# 

# load("models/Arndis/confidence/weighted_confidence_diagnosis_m")
# ps = sample(unique(weighted_confidence_diagnosis_m$data$SubjectPair), 4)
# plot_ppc(weighted_confidence_diagnosis_m,
#           weighted,
#           "posterior",
#           diagnosis = TRUE,
#           random_effects = TRUE,
#           by_participant = TRUE,
#           confidence = TRUE,
#           re_name = "SubjectPair",
#           colour_name = "lOthers",
#           participants = ps)
# 
# load("models/Arndis/confidence/circular_confidence_model")
# ps = sample(unique(circular_confidence_model$data$SubjectPair), 4)
# plot_ppc(circular_confidence_model,
#          circular,
#          "prior",
#          diagnosis = FALSE,
#          random_effects = TRUE,
#          by_participant = TRUE,
#          confidence = TRUE,
#          re_name = "SubjectPair",
#          re_prior = c(.4, 8),
#          colour_name = "lOthers",
#          participants = ps,
#          plot_naive = FALSE)

