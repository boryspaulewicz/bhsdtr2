## -*- coding: utf-8 -*-

## len is the number of columns of the prior matrix
default_prior = function(par, len, prior.par, model, links, K){
    if(prior.par == 'random_scale')
        prior.par = 'fixed_sd'
    unb = unbiased(K)
    Kb2 = round(K / 2)
    priors = list(theta = list(fixed_mu = 0, fixed_sd = log(1.5)),
                  eta = list(fixed_mu = 0, fixed_sd = 3))
    ## prior dla delta zależy od modelu
    if(model == 'metad'){
        priors$delta = list(fixed_mu = rep(acc.to.delta(.75), 2), fixed_sd = rep(.5 * (acc.to.delta(.99) - acc.to.delta(.51)), 2))
    }else{
        priors$delta = list(fixed_mu = acc.to.delta(.75) , fixed_sd = .5 * (acc.to.delta(.99) - acc.to.delta(.51)))
    }
    ## prior dla gamma zależy od funkcji łączącej
    if(links$gamma == 'twoparameter'){
        priors$gamma =  list(fixed_mu = c(0, log(unbiased(K)[K / 2 + 1] - unbiased(K)[K / 2])),
                             fixed_sd = c(priors$eta$fixed_sd, log(2)))
    }else if(links$gamma == 'parsimonious'){
        priors$gamma =  list(fixed_mu = c(0, log(1)),
                             fixed_sd = c(priors$eta$fixed_sd, log(2)))
    }else if(links$gamma == 'softmax'){
        priors$gamma = list(fixed_mu = rep(0, K - 1), fixed_sd = rep(log(100), K - 1))
    }else if(links$gamma == 'log_distance'){
        res.mu = rep(0, K - 1)
        for(k in (Kb2+1):(K - 1))
            res.mu[k] = log(unb[k] - unb[k - 1])
        for(k in (Kb2 - 1):1)
            res.mu[k] = log(unb[k + 1] - unb[k])
        res.sd = rep(log(2), K - 1)
        res.sd[Kb2] = priors$eta$fixed_sd
        priors$gamma = list(fixed_mu = res.mu, fixed_sd = res.sd)
    }else if(links$gamma == 'log_ratio'){
        res.mu = rep(0, K - 1)
        res.sd = rep(log(2), K - 1)
        res.sd[Kb2] = priors$eta$fixed_sd
        res.mu[Kb2 + 1] = log(unb[Kb2 + 1] - unb[Kb2])
        res.mu[Kb2 - 1] = log(1)
        for(k in (Kb2 + 2):(K - 1))
            res.mu[k] = log((unb[k] - unb[k - 1]) / (unb[Kb2 + 1] - unb[Kb2]))
        for(k in (Kb2 - 2):1)
            res.mu[k] = log((unb[k + 1] - unb[k]) / (unb[k + 1] - unb[k]))
        priors$gamma = list(fixed_mu = res.mu, fixed_sd = res.sd)
    }else{
        stop(sprintf('Unknown gamma link %s', links$gamma))
    }
    matrix(priors[[par]][[prior.par]], nrow = par.size(par, model, links, K)[1], ncol = len)
}
