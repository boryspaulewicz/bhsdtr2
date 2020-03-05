## -*- coding: utf-8 -*-

##' Fit a bhsdtr model using jmap or stan
##'
##' @param m a (possibly fitted) bhsdtr model object
##' @param method either 'jmap' (the default) or 'stan'
##' @param stan_optimizations Wether to set the recommended stan
##'     optimizations (mc.cores, auto_write). The default is T.
##' @param jmap_init initial values for the optimizing function (the
##'     default is 0)
##' @param chains The dafault is parallel::detectCores() - 1, i.e.,
##'     one core left for other tasks.
##' @param iter The default is max(20000 / chains, 5000). Only
##'     relevant if method == 'stan'.
##' @param warmup The default is 2000. Only relevant if method ==
##'     'stan'.
##' @param init_r The default is .5. Only relevant if method ==
##'     'stan'.
##' @param ... other arguments passed to the optimizing or stan
##'     function
##' @return a bhsdtr model object with updated $jmap (if method ==
##'     'jmap') or the $stanfit (if method == 'stan') field
##' @export
fit = function(m, method = 'jmap',
               stan_optimizations = T,
               jmap_init = 0,
               chains = parallel::detectCores() - 1, iter = max(20000 / chains, 5000), warmup = 2000, init_r = .5,
               sample.prior = F, ...){
    if(sample.prior)
        method = 'stan'
    if(method == 'jmap'){
        m$jmapfit = optimizing(stan_model(model_code = m$code), m$sdata, init = jmap_init, ...)
        if(m$jmapfit$return_code != 0)
            warning('optimizing did not converge')
    }else if(method == 'stan'){
        stanargs = list(...)
        args = list(chains = chains, iter = iter, init_r = init_r, warmup = warmup)
        for(arg in names(args))
            if(is.null(stanargs[[arg]]))
                stanargs[[arg]] = args[[arg]]
        if(is.null(stanargs$pars)){
            pars = NULL
            for(par in names(m$fixed))
                pars = c(pars, sprintf('%s_fixed', par))
            for(str in c('%s_sd_%d', 'corr_%s_%d', '%s_random_%d'))
                for(par in names(m$random))
                    for(g in 1:length(m$random[[par]]))
                        pars = c(pars, sprintf(str, par, g))
            pars = c(pars, 'counts_new')
            stanargs$pars = m$pars = pars
        }else{
            m$pars = stanargs$pars
        }
        if(stan_optimizations){
            options(mc.cores = parallel::detectCores())
            rstan_options(auto_write = TRUE)
            ## Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
        }
        stanargs$data = m$sdata
        if(sample.prior){
            stanargs$model_code = make.model.code(m$model, m$fixed, m$random, m$links, only_prior = T)
            m$stanfit.prior = do.call(stan, stanargs)
        }else{
            stanargs$model_code = m$code
            m$stanfit = do.call(stan, stanargs)
        }
    }else{
        stop(sprintf('Unknown method %s, must be either jmap or stan', method))
    }
    m
}
