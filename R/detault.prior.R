## -*- coding: utf-8 -*-

## len is the number of columns of the prior matrix, e.g., par == 'delta', prior.par == 'fixed_mu'

##' @export
default.prior = function(m, par, len, prior.par){
    ##! This is not used yet
    mu_of_normal_for_lognormal = function(lognormal_mu, lognormal_sd)log(lognormal_mu^2/sqrt(lognormal_mu^2 + lognormal_sd^2))
    var_of_normal_for_lognormal = function(lognormal_mu, lognormal_sd)log(1 + (lognormal_sd/lognormal_mu)^2)
    model = m$model
    links = m$links
    K = m$sdata$K
    ## Default values for most situations ##! TODO random_scale_lb == ''?
    if(prior.par %in% c('fixed_lb', 'fixed_ub', 'random_scale_lb', 'random_scale_ub'))
        out = ''
    if(prior.par == 'random_nu'){
        out = 1
    }else if(par == 'delta'){
    }else if(par == 'gamma'){
    }else if(par == 'eta'){
    }else if(par == 'theta'){
    }
    prior.par.original = prior.par
    ## Default sd of the prior for the fixed effect is the same as the default sd of the prior for
    ## the random effects' sd
    if(prior.par == 'random_scale')
        prior.par = 'fixed_sd'
    unb = unbiased(K)
    Kb2 = round(K / 2)
    ## Here, we create this list that is later filled with many different priors, even though only
    ## the matrix for the given par and prior.par is returned by this function, because it seems to
    ## be more mangable than the very long switch-case kind of structure that would have to be used
    ## instead
    priors = list(theta = list(fixed_mu = 0, fixed_sd = log(1.5)),
                  eta = list(fixed_mu = 0, fixed_sd = 3),
                  delta = list(), gamma = list())
    priors[[par]][['random_nu']] = 1
    priors[[par]][['fixed_lb']] = priors[[par]][['fixed_ub']] = priors[[par]][['random_scale_lb']] = priors[[par]][['random_scale_ub']] = ''
    if(par == 'theta' & model == 'dpsdt'){
        fixed_mu = -2.2
        if(prior.par.original == 'random_scale'){
            fixed_sd = 4
        }else{
            fixed_sd = 10^6
        }
        priors$theta$fixed_mu = fixed_mu
        priors$theta$fixed_sd = fixed_sd
    }
    ## Lower bounds for id_* link functions. The bounds are always one-dimensional regardless of
    ## the number of dimensions of the parameter.
    if(links[[par]] %in% c('id_log', 'id_ratio_log')){
        if(par == 'delta')
            priors$delta$fixed_lb = 0
        if(par == 'gamma'){
            priors$gamma$fixed_lb = 0
        }
    }
    if('delta' %in% names(links)){
        if(links$delta == 'identity'){
            fixed_mu = exp(acc.to.delta(.75))
            fixed_sd = .5 * (exp(acc.to.delta(.99)) - exp(acc.to.delta(.51)))
        }else if(links$delta == 'id_log'){
            fixed_mu = exp(acc.to.delta(.75))
            if(prior.par.original == 'random_scale'){
                fixed_sd = 4 ## .5 * (acc.to.delta(.99) - acc.to.delta(.51))
            }else{
                fixed_sd = .5 * (exp(acc.to.delta(.99)) - exp(acc.to.delta(.51)))
            }
        }else if(links$delta %in% c('log')){
            fixed_mu = .5 ## acc.to.delta(.75)
            if(prior.par.original == 'random_scale'){
                fixed_sd = 4
            }else{
                fixed_sd = 1 ## .5 * (acc.to.delta(.99) - acc.to.delta(.51))
            }
        }else if(links$delta == 'log_logit'){
            fixed_mu = c(.5, -2.2) ## r ~= .1
            if(prior.par.original == 'random_scale'){
                fixed_sd = 4
            }else{
                fixed_sd = c(1, 10^6) ## .5 * (acc.to.delta(.99) - acc.to.delta(.51))
            }
        }
        priors$delta$fixed_mu = fixed_mu
        priors$delta$fixed_sd = fixed_sd
        ## if(model %in% c('metad'))
        ##     for(prior.par in names(priors$delta))
        ##         priors$delta[[prior.par]] = rep(priors$delta[[prior.par]], 2)
    }
    ## prior for gamma depends on the link function
    if(links$gamma == 'twoparameter'){
        priors$gamma$fixed_mu = c(0, log(unbiased(K)[K / 2 + 1] - unbiased(K)[K / 2]))
        priors$gamma$fixed_sd = c(priors$eta$fixed_sd, log(2))
    }else if(links$gamma == 'parsimonious'){
        priors$gamma$fixed_mu = c(0, log(1))
        priors$gamma$fixed_sd = c(priors$eta$fixed_sd, log(2))
    }else if(links$gamma == 'softmax'){
        priors$gamma$fixed_mu = rep(0, K - 1)
        priors$gamma$fixed_sd = rep(log(100), K - 1)
    }else if(links$gamma == 'log_distance'){
        res.mu = rep(0, K - 1)
        if(K > (Kb2 + 1))
            for(k in (Kb2+1):(K - 1))
                res.mu[k] = log(unb[k] - unb[k - 1])
        if((Kb2 - 1) > 0)
            for(k in (Kb2 - 1):1)
                res.mu[k] = log(unb[k + 1] - unb[k])
        res.sd = rep(log(2), K - 1)
        res.sd[Kb2] = priors$eta$fixed_sd
        priors$gamma$fixed_mu = res.mu
        priors$gamma$fixed_sd = res.sd
    }else if(links$gamma == 'id_log'){
        res.sd = rep(priors$eta$fixed_sd, K - 1)
        if(K > (Kb2 + 1))
            for(k in (Kb2+1):(K - 1))
                res.sd[k] = unb[k] - unb[k - 1]
        if((Kb2 - 1) > 0)
            for(k in (Kb2 - 1):1)
                res.sd[k] = unb[k + 1] - unb[k]
        priors$gamma$fixed_mu = rep(0, K - 1)
        priors$gamma$fixed_sd = res.sd
    }else if(links$gamma == 'log_ratio'){
        res.mu = rep(0, K - 1)
        res.sd = rep(log(2), K - 1)
        res.sd[Kb2] = priors$eta$fixed_sd
        res.mu[Kb2 + 1] = log(unb[Kb2 + 1] - unb[Kb2])
        if(Kb2 > 1)
            res.mu[Kb2 - 1] = log(1)
        if((Kb2 + 2) < K)
            for(k in (Kb2 + 2):(K - 1))
                res.mu[k] = log((unb[k] - unb[k - 1]) / (unb[Kb2 + 1] - unb[Kb2]))
        if((Kb2 - 2) > 0)
            for(k in (Kb2 - 2):1)
                res.mu[k] = log((unb[k + 1] - unb[k]) / (unb[k + 1] - unb[k]))
        priors$gamma$fixed_mu = res.mu
        priors$gamma$fixed_sd = res.sd
    }else{
        stop(sprintf('Unknown gamma link %s', links$gamma))
    }
    ## Bounds are always one-dimensional
    if(prior.par %in% c('fixed_lb', 'fixed_ub')){
        priors[[par]][[prior.par]]
    }else{
        matrix(priors[[par]][[prior.par]], nrow = par.size(par, model, links, K)[1], ncol = len)
    }
}
