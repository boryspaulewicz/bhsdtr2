##' Set the parameters of the prior distributions
##' 
##' @param model a bhsdtr model object
##' @param ... There are two ways of providing the values of the prior
##'     parameters. One is by directly using the names of the prior
##'     parameters, e.g., delta_prior_fixed_mu = log(1.5),
##'     delta_prior_fixed_sd = 2. Note that with the exception of nu
##'     (lkj prior), the prior parameters are internally represented
##'     as matrices of dimension DxC, where D is the dimensionality of
##'     the parameter (delta has 2 in the meta-d' and the dpsdt with
##'     correlated d' and R random effects models, otherwise it has 1,
##'     theta has 1, eta has 1, gamma has K - 1, unless the
##'     twoparameter or the parsimonious link function is used, in
##'     which case it has dim = 2) and C is the number of columns in
##'     the fixed (<par>_prior_fixed_mu, <par>_prior_fixed_sd) or the
##'     random (<par>prior_scale_<g>) effects' model matrix. You can
##'     provide scalars, vectors, or matrices. A vector will fill the
##'     prior parameter matrix in column major order. A more
##'     convenient way is to use arguments of the form: delta =
##'     list(mu = log(2), sd = 4, scale = list('1' = 6), nu = list('1'
##'     = 1000), etc.
##' @return a bhsdtr model object with updated priors
##' @export
set.prior = function(model, ...){
    args = list(...)
    ## set.prior(model, delta = list(...), gamma = list(...))
    if(all(names(args) %in% c('delta', 'gamma', 'theta', 'eta'))){
        ## e.g., delta = list(mu = log(1), scale = list('1' = 3))
        for(par.type in names(args)){
            if(!(par.type %in% c('delta', 'gamma', 'eta', 'theta')))
                stop(sprintf('%s is not a valid parameter name', par.type))
            for(par in names(args[[par.type]])){
                if(!(par %in% c('mu', 'sd', 'lb', 'ub', 'scale', 'nu', 'scale_lb', 'scale_ub')))
                    stop(sprintf('%s %s is not a prior parameter', par.type, par))
                if(par %in% c('scale', 'scale_lb', 'scale_ub')){
                    for(g in names(args[[par.type]][[par]])){
                        par.name = sprintf('%s_prior_%s_%s', par.type, par, g)
                        if(!(par.name %in% names(model$sdata)))
                            stop(sprintf('%s parameter does not exist', par.name))
                        model$sdata[[par.name]] = matrix(args[[par.type]][[par]][[g]],
                                                         nrow = nrow(model$sdata[[sprintf('%s_size', par.type)]]),
                                                         ncol = ncol(model$sdata[[sprintf('Z_%s_ncol_%d', par.type, g)]]))
                    }
                }else if(par == 'nu'){
                    for(g in names(args[[par.type]][[par]])){
                        par.name = sprintf('%s_prior_random_nu_%s', par.type, g)
                        if(!(par.name %in% names(model$sdata)))
                            stop(sprintf('%s parameter does not exist', par.name))
                        model$sdata[[par.name]] = args[[par.type]][[par]][[g]][1]
                    }
                }else{
                    ## fixed effects prior parameters
                    par.name = sprintf('%s_prior_%s', par.type, c(mu = 'fixed_mu', sd = 'fixed_sd')[par])
                    if(!(par.name %in% names(model$sdata)))
                        stop(sprintf('%s parameter does not exist', par.name))
                    model$sdata[[par.name]] = matrix(args[[par.type]][[par]],
                                                     nrow = nrow(model$sdata[[par.name]]),
                                                     ncol = ncol(model$sdata[[par.name]]))
                }
            }
        }
    }else{
        for(arg in names(args)){
            if(length(grep('prior', arg)) == 0)
                stop(sprintf('%s is not a prior parameter', arg))
            if(is.null(model$sdata[[arg]])){
                stop(sprintf('Prior parameter %s does not exist in this model', arg))
            }else{
                model$sdata[[arg]] = matrix(args[[arg]], nrow = nrow(model$sdata[[arg]]), ncol = ncol(model$sdata[[arg]]))
            }
        }
    }
    model$code = parse.model.code(model)
    model
}
