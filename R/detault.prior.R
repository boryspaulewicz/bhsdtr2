## -*- coding: utf-8 -*-

## len is the number of columns of the prior matrix
default_prior = function(par, len, prior.par, model, links, K){
    if(prior.par == 'random_scale')
        prior.par = 'fixed_sd'
    priors = list(theta = list(fixed_mu = 0, fixed_sd = 1),
                  eta = list(fixed_mu = 0, fixed_sd = 10))
    ## prior dla delta zależy od modelu
    if(model == 'metad'){
        priors$delta = list(fixed_mu = c(0, 0), fixed_sd = c(2, 2))
    }else{
        priors$delta = list(fixed_mu = 0, fixed_sd = 2)
    }
    ## prior dla gamma zależy od funkcji łączącej
    if(links$gamma %in% c('twoparameter', 'parsimonious')){
        priors$gamma =  list(fixed_mu = c(0, log(2 / K)),
                             fixed_sd = c(10, log(2)))
    }else if(links$gamma == 'softmax'){
        priors$gamma = list(fixed_mu = rep(0, K - 1), fixed_sd = rep(log(100), K - 1))
    }else{
        ## Uproszczenie - trzeba inaczej modelować crit, log(dist) i log(ratio(dist))
        priors$gamma = list(fixed_mu = rep(0, K - 1), fixed_sd = rep(log(2), K - 1))
    }
    matrix(priors[[par]][[prior.par]], nrow = par.size(par, model, links, K)[1], ncol = len)
}
