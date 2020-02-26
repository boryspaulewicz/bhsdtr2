##' ## -*- coding: utf-8 -*-

##' fit a bhsdtr model using ml or stan
##' 
##' @export
fit = function(m, method = 'ml',
               ml_init = 0,
               chains = parallel::detectCores() - 1, iter = max(20000 / chains, 5000), warmup = 2000, init_r = .5,
               stan_optimizations = T,
               ...){
    if(method == 'ml'){
        m$mlfit = optimizing(stan_model(model_code = m$code), m$sdata, init = ml_init, ...)
        if(m$mlfit$return_code != 0)
            warning('optimizing did not converge')
    }
    if(method == 'stan'){
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
            Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
        }
        stanargs$model_code = m$code
        stanargs$data = m$sdata
        m$stanfit = do.call(stan, stanargs)
    }
    m
}
