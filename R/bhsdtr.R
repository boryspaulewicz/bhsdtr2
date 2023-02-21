## -*- coding: utf-8 -*-

##' Creating and fitting ordinal models with ordered thresholds
##'
##' This is the main method of fitting ordinal models in the bhsdtr2
##' package. The model type (e.g., SDT, UV SDT, DPSDT, meta-d') is
##' inferred automatically based on the supplied model formulae. An EV
##' SDT model is fitted if the user provides model formulae for dprim
##' and thr, a meta-d' model is fitted if there is a formula for the
##' metad parameter vector, an UV version is fitted if there is
##' formula for the sd_ratio parameter. Fitting an DPSDT model is a
##' bit more complicated because there are two different versions of
##' this model and because the probability of recall may be non-zero
##' only for one stimulus class (e.g., "single-interval" tasks) or for
##' both stimulus classes (2AFC tasks). In order to fit an DPSDT model
##' the user has to introduce a recall indicator in the
##' response_formula, e.g., resp ~ stim + possible_recall, where
##' possible_recall is 1 iff recall is possible (a correct response
##' with maximum confidence) on that trial. If d' and R may be
##' correlated, then the user may model the dprimr parameter vector
##' (the first element is the d' parameter, the second element is R),
##' otherwise three formulas have to be provided: one for dprim, one
##' for thr, and one for R. At present, only the logit link function
##' is implemented for the R parameter.
##' 
##' @param model_formulae a vector of two or more lmer-style model formulae for the parameters of
##'     the model. Every such vector must contain a model formula for the thr (=thresholds)
##'     parameter, and a model formula for either the dprim, the metad, or the mean (general ordinal
##'     model) parameter. This vector may have three elements iff it also contains a formula for the
##'     sdratio (Unequal Variance models) parameter.
##' @param response_formula a formula specifying the names of the response variable and the stimulus
##'     variable. If this is an SDT-like model this must be of the form resp ~ stim [+ modifier],
##'     where resp is the name of the combined response variable or the binary response variable,
##'     stim is the name of the (binary) stimulus variable, and the optional modifier indicates
##'     something model-specific (e.g., the possibility of recall in DPSDT).
##' @param links a list specyfing the link functions, the default is list(gamma =
##'     'log_distance'). This list will be filled with the default link function specifications for
##'     the remaining parameters if these are not specified. Currently you can set gamma to
##'     'log_distance' (the default), 'log_ratio', 'twoparameter', 'parsimonious', and
##'     'softmax'. The log_distance and log_ratio link functions are described in the Readme file in
##'     the github repository, the parsimonious link function represents the thresholds using two
##'     parameters, the unconstrained main criterion and the logarithm of the scaling factor which
##'     is used to stretch the "unbiased" thresholds so that the distribution of responses fits. The
##'     softmax link function is described in the bhsdtr preprint. You can also choose between delta
##'     = 'log' (the default, which covers both the d' and the meta-d' parameters), and delta =
##'     'identity' (i.e., no conversion between delta and d' or meta-d').
##' @param method [= 'jmap'] a string which specifies how to fit the model. The default is 'jmap',
##'     which results in an attempt to quickly maximize the joint posterior by using the
##'     rstan::optimizing function, but you can also choose 'stan', or, if you do not want the model
##'     to be fitted right away, you can use anything else.
##' @param prior [= NULL] a list with non-default prior settings which will be passed to the
##'     set.prior function. See the set.prior function documentation for details.
##' @param thresholds_scale thresholds scaling factor for the softmax gamma link function (see the
##'     bhsdtr preprint).
##' @param force.id_log If FALSE (the default) separate intercepts parametrization will be required
##'     for the delta / d' / meta-d' and the gamma / thresholds' model matrices when using the
##'     id_log link function.
##' @param fit If TRUE (the default) then also fit the model
##' @param ... arguments to be passed to the stan function.
##' @return a bhsdtr_model object, which is an S3 class object. If it was fitted using the jmap
##'     (stan) method, $jmapfit ($stanfit) contains the result. You can refit the model using either
##'     method (fit(model, ...)) and the $jmapfit or the $stanfit slot will be updated.
##' @examples
##' \donttest{
##'     gabor$r = combined.response(gabor$stim, gabor$rating, gabor$acc)
##'     ## A hierarchical SDT model
##'     m.jmap = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id)),
##'                   r ~ stim,
##'                   gabor)
##' 
##'     ## Posterior samples (here single samples = maximum joint
##'     ## posterior point estimates) of d' for every unique combination
##'     ## of the predictors, here duration and order, because
##'     ## dprim ~ duration * order.
##'     samples(m.jmap, 'dprim')
##' 
##'     ## Thresholds
##'     samples(m.jmap, 'thr')
##' 
##'     ## Fitting the model using the Stan sampler
##'     m.stan = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id)),
##'                     r ~ stim,
##'                     gabor, method = 'stan')
##'     (smp = samples(m.stan, 'thr'))
##' 
##'     ## This is how you can access the stanfit object
##'     print(m.stan$stanfit, probs = c(.025, .975), pars = c('delta_fixed', 'gamma_fixed'))
##' 
##'     ## The plot method lets you see if the model fits. Here we
##'     ## specify the variables (duration and order) to see only
##'     ## the duration x order panels.
##'     library(ggplot2)
##'     plot(m.stan, vs = c('duration', 'order'), verbose = F)
##' 
##'     ## A simple contrast calculated on the posterior threshold samples
##'     round(t(apply(smp, 2, function(x)quantile(x[,'DECISION-RATING'] - x[,'RATING-DECISION'],
##'                                               c(.025, .975)))), 2)
##' }
##' @export
bhsdtr = function(model_formulae, response_formula, data,
                  links = list(gamma = 'log_distance'), method = 'jmap',
                  prior = list(), sample.prior = F, thresholds_scale = 2, force.id_log = F, fit = T, ...){
    if(method == 'ml'){
        method = 'jmap'
        stop('method = \'ml\' is no longer a valid argument, use method = \'jmap\' instead.
The fitted object will be stored in the $jmapfit field of the bhsdtr model object')
    }
    data = as.data.frame(data)
    ## The model object is filled with parsed model formulae and original parameter names (e.g., 'dprim')
    m = parse.model_formulae(model_formulae, data)
    ## We can now identify the variables that have to preserved when aggregating the data
    vnames = NULL
    for(p in names(m$fixed))
        vnames = c(vnames, names(get_all_vars(m$fixed[[p]], data)))
    for(p in names(m$random))
        for(g in 1:length(m$random[[p]]))
            for(k in c('model.formula', 'group.formula'))
                vnames = c(vnames, names(get_all_vars(m$random[[p]][[g]][[k]], data)))
    vnames = unique(vnames)
    ## Inferring the model type
    if(setequal(m$pars, c('dprim', 'thr', 'sdratio'))){
        m$model = 'uvsdt'
    }else if(setequal(m$pars, c('dprim', 'thr'))){
        m$model = 'sdt'
    }else if(setequal(m$pars, c('metad', 'thr'))){
        m$model = 'metad'
    }else if(setequal(m$pars, c('metad', 'thr', 'sdratio'))){
        m$model = 'uvmetad'
    }else if(setequal(m$pars, c('dprimr', 'thr'))){
        m$model = 'dpsdtcor'
    }else if(setequal(m$pars, c('dprim', 'thr', 'R'))){
        m$model = 'dpsdt'
    }else if(setequal(m$pars, c('dprim', 'thr', 'r'))){
        m$model = 'dpsdt'
    }else if(setequal(m$pars, c('mean', 'thr'))){
        m$model = 'ordinal'
    }else{
        stop(sprintf('Unknown model. Parameters: %s', paste(pars, collapse = ' ')))
    }
    ## Creating resp and stim variables: only resp ~ stim [+ modifier] (sdt family)
    ## or resp ~ 1 (ordinal family) allowed
    resp.stim.model.frame = model.frame(response_formula, data)
    resp.stim.vnames = names(get_all_vars(response_formula, data))
    if(ncol(resp.stim.model.frame) > 3)
        stop(sprintf('Response formula must be of the form response_variable ~ stimulus_variable [+ modifier] or response_variable ~ 1: %s',
                     as.character(response_formula)))
    ## We do not allow for gaps in the response variable e.g., 1, 1,
    ## 2, 4, 5 (no 3s) nor in the stimulus variable: all integers
    ## between 1 and max(stim) or max(resp) have to be present.
    resp.stim.vars = list(resp = fix.index.gaps(resp.stim.model.frame[, 1], names(resp.stim.vnames)[1], T))
    if(length(resp.stim.vnames) > 1){
        resp.stim.vars$stim = fix.index.gaps(resp.stim.model.frame[, 2], resp.stim.vnames[2], T)
        if(length(resp.stim.vnames) > 2)
            resp.stim.vars$modifier = resp.stim.model.frame[, 3]
    }else{
        ## dummy stimulus variable (response_formula was of the form
        ## resp ~ 1)
        resp.stim.vars$stim = rep(1, nrow(resp.stim.model.frame))
    }
    K = max(resp.stim.vars$resp, na.rm = T)
    ## Original dataset
    m$data = data
    ## Agreggation: only the relevant variables without the response
    ## and the stimulus variables. If there is only one variable, we
    ## have to make sure data does not become a vector.
    data = data[, vnames, drop = F]
    res = aggregate.data(data, resp.stim.vars, K)
    m$adata = res$adata
    ## checking if the link functions are valid, inferring the default
    ## link functions for omitted parameters
    m$links = fill.model.links(m, links, adata$data, K)
    ## Creating data structures and model code for stan
    m$sdata = list(N = nrow(m$adata$data),
                              K = K,
                              Kb2 = round(K / 2),
                              PRINT = 0,
                              unbiased = fix.stan.dim(unbiased(K)),
                              thresholds_scale = thresholds_scale,
                              counts = res$counts,
                              modifier = m$adata$modifier)
    rm(res)
    ## in SDT models stim_sign = -1, 1 (this variable slightly simplifies the model code)
    m$sdata$stim_sign = 2 * as.numeric(as.factor(as.character(m$adata$stimulus))) - 3
    m$sdata = make.sdata.matrices(m, force.id_log)
    if(m$model == 'ordinal')
        m$sdata$eta_is_fixed[1,1] = 1
    ## bhsdtr_model object
    m = append(m, list(resp = resp.stim.vars$resp, stim = resp.stim.vars$stim))
    class(m) = c('bhsdtr_model', class(m))
    if(!is.null(prior)){
        prior$model = m
        m = do.call(set.prior, prior)
    }
    m$code = parse.model.code(m)
    if(fit){
        if(sample.prior){
            method = 'stan'
            m = fit(m, method, ...)
            m = fit(m, method, sample.prior = T, ...)
        }else if(method %in% c('jmap', 'stan'))
            m = fit(m, method, ...)
    }
    m
}

parse.model_formulae = function(model_formulae, data){
    fixed = random = list()
    pars = NULL
    if(length(model_formulae) < 2)
        stop('Vector of model formulae must be of length > 1')
    for(i in 1:length(model_formulae)){
        ## fixed effects formula
        fixed.formula = lme4::nobars(model_formulae[[i]])
        ## delta, gamma, ...
        pars = c(pars, as.character(fixed.formula[[2]]))
        ## dprim -> delta, etc
        par = par.to.linked(as.character(fixed.formula[[2]]))
        ## ~ <fixed effects part> ...
        fixed[[par]] = lme4::nobars(fixed.formula)[-2]
        ## random effects part, e.g., (1 | id)
        ranef.formulae = lme4::findbars(model_formulae[[i]])
        if(!is.null(ranef.formulae)){
            random[[par]] = list()
            for(g in 1:length(ranef.formulae)){
                ## ~ (duration | id) -> ~ duration
                random[[par]][[g]] = c(model.formula = as.formula(paste('~', as.character(ranef.formulae[[g]])[2])))
                ## ~ (duration | id) -> ~ id
                random[[par]][[g]]$group.formula = as.formula(paste('~', as.character(ranef.formulae[[g]])[3]))
                ranef.model.frame = model.frame(random[[par]][[g]]$group.formula, data)
                ## max(id)
                random[[par]][[g]]$group.size = max(fix.index.gaps(ranef.model.frame[, 1], names(ranef.model.frame)[1], T))
                ## 'id'
                random[[par]][[g]]$group.name = names(get_all_vars(random[[par]][[g]]$group.formula, data))
                ## just checking
                if(length(random[[par]][[g]]$group.name) > 1)
                    stop(sprintf('More than one grouping variable in %s %s | %s', par, as.character(random[[par]][[g]]$model.formula),
                                 paste(random[[par]][[g]]$group.name, collapse = ' ')))
            }
        }
    }
    list(fixed = fixed, random = random, pars = pars)
}

## data structures required by the model
make.sdata.matrices = function(m, force.id_log = F){
    ## Fixed effects model matrices and priors
    for(par in names(m$fixed)){
        v = sprintf('X_%s', par)
        m$sdata[[v]] = model.matrix(m$fixed[[par]], m$adata$data)
        if(m$links[[par]] == 'id_log')
            if(!is.separate.intercepts(m$sdata[[v]]) & !force.id_log)
                stop(sprintf('id_log link requires separate intercepts parametrization (e.g., ~ -1 + f1:f2),
  found: %s %s\n',
                             par, paste(as.character(m$fixed[[par]]), collapse = ' ')))
        m$sdata[[sprintf('X_%s_ncol', par)]] = ncol(m$sdata[[v]])
        m$sdata[[sprintf('%s_is_fixed', par)]] = m$sdata[[sprintf('%s_fixed_value', par)]] =
            matrix(0, nrow = par.size(par, m$model, m$links, m$sdata$K)[1], ncol = ncol(m$sdata[[v]]))
        ## Fixed effects priors
        for(prior.bound in c('fixed_lb', 'fixed_ub'))
            m$sdata[[sprintf('%s_prior_%s', par, prior.bound)]] =
                default.prior(m, par, ncol(m$sdata[[v]]), prior.bound)
        for(prior.par in c('fixed_mu', 'fixed_sd'))
            m$sdata[[sprintf('%s_prior_%s', par, prior.par)]] = default.prior(m, par, ncol(m$sdata[[v]]), prior.par)
    }
    ## Random effects model matrices and priors
    for(par in names(m$random)){
        Z = g = g.original = list()
        Z_ncol = g_max = NULL
        for(i in 1:length(m$random[[par]])){
            Z[[i]] = model.matrix(m$random[[par]][[i]]$model.formula, m$adata$data)
            if(m$links[[par]] == 'id_log')
                if(!is.separate.intercepts(Z[[i]]) & !force.id_log)
                    stop(sprintf('id_log link requires separate intercepts parametrization (e.g., ~ -1 + f1:f2),
  found: %s %s | %s',
  par, paste(as.character(m$random[[par]][[i]]$model.formula), collapse = ' '),
  m$random[[par]][[i]]$group.name))
            Z_ncol = c(Z_ncol, ncol(Z[[i]]))
            ## group indicator variable
            mf = model.frame(m$random[[par]][[i]]$group.formula, m$adata$data)
            if(ncol(mf) == 1){
                g.original[[i]] = mf[, 1]
                g[[i]] = fix.index.gaps(mf[, 1], names(mf)[1], T)
                g_max = c(g_max, max(g[[i]]))
            }else{
                stop(sprintf('Something is wrong with the grouping factor for %s: %s',
                             par, paste(as.character(random[[par]][[i]]$group.formula, collapse = ' '))))
            }
        }
        ## Because random effects' matrices may differ in the number
        ## of columns we have to use indexed names, e.g., Z_delta_1,
        ## Z_delta_2, etc.
        for(i in 1:length(Z)){
            m$sdata[[sprintf('%s_group_max_%d', par, i)]] = g_max[i]
            m$sdata[[sprintf('%s_group_%d', par, i)]] = g[[i]]
            m$sdata[[sprintf('%s_group_%d_original', par, i)]] = g.original[[i]]
            m$sdata[[sprintf('Z_%s_%d', par, i)]] = Z[[i]]
            m$sdata[[sprintf('Z_%s_ncol_%d', par, i)]] = Z_ncol[i]
            m$sdata[[sprintf('zeros_%s_%d', par, i)]] = fix.stan.dim(rep(0, par.size(par, m$model, m$links, m$sdata$K)[1] * Z_ncol[i]))
            ## scale and nu priors
            m$sdata[[sprintf('%s_prior_random_nu_%d', par, i)]] = 1
            m$sdata[[sprintf('%s_prior_random_sd_lb_%d', par, i)]] = 0
            m$sdata[[sprintf('%s_prior_random_sd_ub_%d', par, i)]] = NA
            m$sdata[[sprintf('%s_prior_random_scale_%d', par, i)]] = default.prior(m, par, Z_ncol[i], 'random_scale')
        }
    }
    ## Parameter matrix dimensions for par and par_ = link(par)
    for(par in unique(c(names(m$fixed), names(m$random)))){
        m$sdata[[sprintf('%s_size', par)]] = par.size(par, m$model, m$links, m$sdata$K)[1]
        m$sdata[[sprintf('%s_size_', par)]] = par.size(par, m$model, m$links, m$sdata$K)[2]
    }
    m$sdata
}

aggregate.data = function(data, resp.stim.vars, K){
    vnames = names(data)
    ## Unique stimulus and response names for aggregation
    resp.stim.unique.names = c(resp = NA, stim = NA, modifier = NULL)
    for(i in 1:length(resp.stim.vars)){
        ## unique names for resp and stim, _ is safe, some other
        ## choices (e.g., @) result in error
        resp.stim.unique.names[names(resp.stim.vars)[i]] = paste(c(names(resp.stim.vars)[i], vnames), collapse = '_')
        ## add uniquely named resp and stim and modifier vars to data
        data[[resp.stim.unique.names[i]]] = resp.stim.vars[[i]]
    }
    ## Removing NA rows
    for(v in names(data))
        data = data[!is.na(data[[v]]),]
    ## aggregation
    res = plyr::ddply(data, unique(c(vnames, resp.stim.unique.names[names(resp.stim.vars) == 'stim'],
                                     resp.stim.unique.names[names(resp.stim.vars) == 'modifier'])),
                      ## We are adding the 1:K vector to make sure
                      ## that the counts for every k are here, the -1
                      ## term corrects for this
                      function(df)table(c(df[[resp.stim.unique.names[names(resp.stim.vars) == 'resp']]], 1:K)) - 1)
    counts = res[, c((ncol(res) - K + 1):ncol(res))]
    ## aggregated data object
    adata = list(data = res[, setdiff(vnames, resp.stim.unique.names), drop = F])
    adata$stimulus = res[[resp.stim.unique.names[names(resp.stim.vars) == 'stim']]]
    if(length(resp.stim.vars) == 3){
        adata$modifier = res[[resp.stim.unique.names[names(resp.stim.vars) == 'modifier']]]
    }else{
        adata$modifier = rep(1, nrow(adata$data))
    }
    list(adata = adata, counts = counts)
}

fill.model.links = function(model, links, data, K){
    model_pars = list(sdt = c('delta', 'gamma'),
                      metad = c('delta', 'gamma'),
                      dpsdtcor = c('delta', 'gamma'),
                      dpsdt = c('delta', 'gamma', 'theta'),
                      uvmetad = c('delta', 'gamma', 'theta'),
                      uvsdt = c('delta', 'gamma', 'theta'),
                      ordinal = c('eta', 'gamma'))
                      ## uvordinal = c('eta', 'gamma', 'theta'))
    pars = unique(c(names(model$fixed), names(model$random)))
    if(!setequal(model_pars[[model$model]], pars))
        stop(sprintf('Model %s must specify effects for: %s', model, paste(model_pars[[model$model]], collapse = ', ')))
    if(!all(names(links) %in% pars))
        stop(sprintf('Found %s link specification but no model formula(e)',
                     paste(names(links)[!(names(links) %in% pars)], collapse = ', ')))
    ## data df = 2K - 2, uvsdt df = K - 1 + 2, 2K - 2 = K + 1 -> K = 3
    if((K < 3) & (model$model %in% c('uvsdt', 'metad', 'dpsdtcor')))
        stop(sprintf('Model %s needs K > 2', model))
    ## The first link function in the vector is the default one
    par_links = list(gamma = c('log_distance', 'log_ratio', 'softmax', 'twoparameter', 'parsimonious', 'identity', 'id_log'),
                     delta = c('log', 'id_ratio_log', 'id_log', 'log_logit', 'identity'),
                     theta = c('log', 'logit'),
                     eta = c('identity'))
    for(par in names(links))
       if(!(links[[par]] %in% par_links[[par]]))
            stop(sprintf('Link %s not allowed for %s', links[[par]], par))
    ## Fill the missing fields in the links list
    for(par in model_pars[[model$model]])
        if(!(par %in% names(links))){
            links[[par]] = par_links[[par]][1]
        }
    ## For some models some links are forced
    if(model$model == 'dpsdt')
        links$theta = 'logit'
    if(model$model == 'dpsdtcor')
        links$delta = 'log_logit'
    links
}
## Ok

## parameter dimensionality (pre-link and post-link)
par.size = function(par, model, links, K){
    if(par == 'gamma'){
        if(links$gamma %in% c('twoparameter', 'parsimonious')){
            ## gamma_size = 2, gamma_size_ = K - 1
            s = c(2, K - 1)
        }else{
            s = rep(K - 1, 2)
        }
    }else{
        s = list(sdt = c(delta = 1),
                 uvsdt = c(delta = 1, theta = 1),
                 metad = c(delta = 2),
                 dpsdtcor = c(delta = 2),
                 dpsdt = c(delta = 1, theta = 1),
                 uvmetad = c(delta = 2, theta = 1),
                 ordinal = c(eta = 1),
                 uvordinal = c(eta = 1, theta = 1))[[model]][par]
        s = rep(s, 2)
    }
    if(any(is.na(s)))
        stop(sprintf('Cannot determine parameter dim: par = %s, model = %s, link = %s, K = %d', par, model, links, K))
    s
}
## Ok
