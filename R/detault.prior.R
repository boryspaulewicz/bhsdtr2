## -*- coding: utf-8 -*-

## len is the number of columns of the prior matrix
default.prior = function(par, len, prior.par, model, links, K){
    prior.par.original = prior.par
    if(prior.par == 'random_scale')
        prior.par = 'fixed_sd'
    unb = unbiased(K)
    Kb2 = round(K / 2)
    priors = list(theta = list(fixed_mu = 0, fixed_sd = log(1.5)),
                  eta = list(fixed_mu = 0, fixed_sd = 3))
    ## prior dla delta zależy od modelu i funkcji łączącej
    if(links$delta == 'identity'){
        fixed_mu = exp(acc.to.delta(.75))
        fixed_sd = .5 * (exp(acc.to.delta(.99)) - exp(acc.to.delta(.51)))
    }else if(links$delta == 'id_log'){
        fixed_mu = exp(acc.to.delta(.75))
        if(prior.par.original == 'random_scale'){
            fixed_sd = 6 ## .5 * (acc.to.delta(.99) - acc.to.delta(.51))
        }else{
            fixed_sd = .5 * (exp(acc.to.delta(.99)) - exp(acc.to.delta(.51)))
        }
    }else if(links$delta == 'log'){
        fixed_mu = .5 ## acc.to.delta(.75)
        if(prior.par.original == 'random_scale'){
            fixed_sd = 6
        }else{
            fixed_sd = 1 ## .5 * (acc.to.delta(.99) - acc.to.delta(.51))
        }
    }
    if(model == 'metad'){
        size = 2
    }else{
        size = 1
    }
    priors$delta = list(fixed_mu = rep(fixed_mu, size), fixed_sd = rep(fixed_sd, size))
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
        if(K > (Kb2 + 1))
            for(k in (Kb2+1):(K - 1))
                res.mu[k] = log(unb[k] - unb[k - 1])
        if((Kb2 - 1) > 0)
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
        if(Kb2 > 1)
            res.mu[Kb2 - 1] = log(1)
        if((Kb2 + 2) < K)
            for(k in (Kb2 + 2):(K - 1))
                res.mu[k] = log((unb[k] - unb[k - 1]) / (unb[Kb2 + 1] - unb[Kb2]))
        if((Kb2 - 2) > 0)
            for(k in (Kb2 - 2):1)
                res.mu[k] = log((unb[k + 1] - unb[k]) / (unb[k + 1] - unb[k]))
        priors$gamma = list(fixed_mu = res.mu, fixed_sd = res.sd)
    }else{
        stop(sprintf('Unknown gamma link %s', links$gamma))
    }
    matrix(priors[[par]][[prior.par]], nrow = par.size(par, model, links, K)[1], ncol = len)
}
