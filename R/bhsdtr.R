## -*- coding: utf-8 -*-

##' main interface
##' 
##' @export
bhsdtr = function(model_formulae, response_formula, data,
                  links = list(gamma = 'log_distance'), fit_method = 'ml',
                  thresholds_scale = 2, ...){
    fixed = random = list()
    data = as.data.frame(data)
    vnames = pars = NULL
    if(length(model_formulae) < 2)
        stop('Vector of model formulae must be of length > 1')
    for(i in 1:length(model_formulae)){
        ## fixed effects formula
        fixed.formula = lme4::nobars(model_formulae[[i]])
        pars = c(pars, as.character(fixed.formula[[2]]))
        ## dprim -> delta, etc
        par = par.to.linked(as.character(fixed.formula[[2]]))
        ## ~ <fixed effects part> ...
        fixed[[par]] = lme4::nobars(fixed.formula)[-2]
        ## storing variable names for aggregation
        vnames = c(vnames, names(get_all_vars(fixed.formula[-2], data)))
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
                ## updating vector of variable names for aggregation
                for(k in c('model.formula', 'group.formula'))
                    vnames = c(vnames, names(get_all_vars(random[[par]][[g]][[k]], data)))
            }
        }
    }
    vnames = unique(vnames)
    if(setequal(pars, c('dprim', 'thr', 'sdratio'))){
        model = 'uvsdt'
    }else if(setequal(pars, c('dprim', 'thr'))){
        model = 'sdt'
    }else if(setequal(pars, c('metad', 'thr'))){
        model = 'metad'
    }else if(setequal(pars, c('mean', 'thr'))){
        model = 'ordinal'
    }else{
        stop(sprintf('Unknown model. Parameters: %s', paste(pars, collapse = ' ')))
    }
    ## Creating resp and stim variables: only resp ~ stim (sdt family)
    ## or resp ~ 1 (ordinal family) allowed
    resp.stim.model.frame = model.frame(response_formula, data)
    resp.stim.vnames = names(get_all_vars(response_formula, data))
    if(ncol(resp.stim.model.frame) > 2)
        stop(sprintf('Response formula must be of the form response_variable ~ stimulus_variable or response_variable ~ 1: %s',
                     as.character(response_formula)))
    ## We do not allow for gaps in the response variable e.g., 1, 1,
    ## 2, 4, 5 (no 3s) nor in the stimulus variable: all integers
    ## between 1 and max(stim) or max(resp) have to be present.
    resp.stim.vars = list(resp = fix.index.gaps(resp.stim.model.frame[, 1], names(resp.stim.vnames)[1], T))
    if(length(resp.stim.vnames) == 2){
        resp.stim.vars$stim = fix.index.gaps(resp.stim.model.frame[, 2], resp.stim.vnames[2], T)
    }else{
        ## dummy stimulus variable (response_formula was of the form
        ## resp ~ 1)
        resp.stim.vars$stim = rep(1, nrow(resp.stim.model.frame))
    }
    ## Agreggation
    K = max(resp.stim.vars$resp, na.rm = T)
    ## Only the relevant variables without response and stimulus variables. 
    ## If there is only one variable, we have to make sure data does not
    ## become a vector
    data = data[, vnames, drop = F]
    res = aggregate.data(data, resp.stim.vars, K)
    adata = res$adata
    links = fix.model.links(model, fixed, random, links, K)
    ## Creatinng data structures and model code for stan
    sdata = list(N = nrow(adata$data), K = K, Kb2 = round(K / 2), PRINT = 0,
                 unbiased = fix.stan.dim(unbiased(K)), thresholds_scale = thresholds_scale,
                 counts = res$counts)
    rm(res)
    ## in SDT models stim_sign = -1, 1 (this variable slightly simplifies the model code)
    sdata$stim_sign = 2 * as.numeric(as.factor(as.character(adata$stimulus))) - 3
    sdata = sdata.matrices(sdata, adata, fixed, random, model, links)
    if(model == 'ordinal')
        sdata$eta_is_fixed[1,1] = 1
    ## Model code
    par_types = unique(names(fixed), names(random))
    parsed = parse.PAR(readLines(stan.file('ordinal_template.stan')), par_types)
    parsed = parse.likelihood(parsed, model)
    parsed = parse.random(parsed, random, par_types)
    for(par in names(links))
        parsed = parse.link(parsed, par, links[[par]])
    ## bhsdtr object
    m = list(fixed = fixed, random = random, adata = adata, sdata = sdata, model = model, links = links, data_size = nrow(data),
             code = paste(parsed, collapse = '\n'))
    class(m) = c('bhsdtr_model', class(m))
    if(fit_method %in% c('ml', 'stan'))
        m = fit(m, fit_method, ...)
    m
}

sdata.matrices = function(sdata, adata, fixed, random, model, links){
    ## Fixed effects' model matrices
    for(par in names(fixed)){
        v = sprintf('X_%s', par)
        sdata[[v]] = model.matrix(fixed[[par]], adata$data)
        sdata[[sprintf('X_%s_ncol', par)]] = ncol(sdata[[v]])
        sdata[[sprintf('%s_is_fixed', par)]] = sdata[[sprintf('%s_fixed_value', par)]] =
            matrix(0, nrow = par.size(par, model, links, sdata$K)[1], ncol = ncol(sdata[[v]]))
        ## priors
        for(prior.par in c('fixed_mu', 'fixed_sd'))
            sdata[[sprintf('%s_prior_%s', par, prior.par)]] = default_prior(par, ncol(sdata[[v]]), prior.par, model, links, sdata$K)
    }
    ## Random effects' model matrices
    for(par in names(random)){
        Z = g = list()
        Z_ncol = g_max = NULL
        for(i in 1:length(random[[par]])){
            Z[[i]] = model.matrix(random[[par]][[i]]$model.formula, adata$data)
            Z_ncol = c(Z_ncol, ncol(Z[[i]]))
            ## group indicator variable
            mf = model.frame(random[[par]][[i]]$group.formula, adata$data)
            if(ncol(mf) == 1){
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
            sdata[[sprintf('%s_group_max_%d', par, i)]] = g_max[i]
            sdata[[sprintf('%s_group_%d', par, i)]] = g[[i]]
            sdata[[sprintf('Z_%s_%d', par, i)]] = Z[[i]]
            sdata[[sprintf('Z_%s_ncol_%d', par, i)]] = Z_ncol[i]
            ## scale and nu priors
            sdata[[sprintf('%s_prior_nu_%d', par, i)]] = 1
            sdata[[sprintf('%s_prior_scale_%d', par, i)]] = default_prior(par, Z_ncol[i], 'random_scale', model, links, sdata$K)
        }
    }
    ## Parameter matrix dimensions for par and par_ = link(par)
    for(par in unique(c(names(fixed), names(random)))){
        sdata[[sprintf('%s_size', par)]] = par.size(par, model, links, sdata$K)[1]
        sdata[[sprintf('%s_size_', par)]] = par.size(par, model, links, sdata$K)[2]
    }
    sdata
}

aggregate.data = function(data, resp.stim.vars, K){
    vnames = names(data)
    ## Unique stimulus and response names for aggregation
    resp.stim.unique.names = c(resp = NA, stim = NA)
    for(i in 1:2){
        ## unique names for resp and stim, _ is safe, some other
        ## choices (e.g., @) result in error
        resp.stim.unique.names[names(resp.stim.vars)[i]] = paste(c(names(resp.stim.vars)[i], vnames), collapse = '_')
        ## add uniquely named resp and stim vars to data
        data[[resp.stim.unique.names[i]]] = resp.stim.vars[[i]]
    }
    ## aggregation
    res = plyr::ddply(data, unique(c(vnames, resp.stim.unique.names[names(resp.stim.vars) == 'stim'])),
                      ## We are adding the 1:K vector to make sure
                      ## that the counts for every k are here, the -1
                      ## term corrects for this
                      function(df)table(c(df[[resp.stim.unique.names[names(resp.stim.vars) == 'resp']]], 1:K)) - 1)
    counts = res[, c((ncol(res) - K + 1):ncol(res))]
    ## aggregated data object
    adata = list(data = res[, setdiff(vnames, resp.stim.unique.names), drop = F])
    adata$stimulus = res[[resp.stim.unique.names[names(resp.stim.vars) == 'stim']]]
    list(adata = adata, counts = counts)
}

fix.model.links = function(model, fixed, random, links, K){
    model_pars = list(sdt = c('delta', 'gamma'),
                      metad = c('delta', 'gamma'),
                      uvsdt = c('delta', 'gamma', 'theta'),
                      ordinal = c('eta', 'gamma'))
                      ## uvordinal = c('eta', 'gamma', 'theta'))
    pars = unique(c(names(fixed), names(random)))
    if(!setequal(model_pars[[model]], pars))
        stop(sprintf('Model %s must specify effects for: %s', model, paste(model_pars[[model]], collapse = ', ')))
    if(!(names(links) %in% pars))
        stop(sprintf('Found %s link specification but no model formula(e)',
                     paste(names(links)[!(names(links) %in% pars)], collapse = ', ')))
    ## data df = 2K - 2, uvsdt df = K - 1 + 2, 2K - 2 = K + 1 -> K = 3
    if((K < 3) & (model %in% c('uvsdt', 'metad')))
        stop(sprintf('Model %s needs K > 2', model))
    par_links = list(gamma = c('log_distance', 'log_ratio', 'softmax', 'twoparameter', 'parsimonious', 'identity'),
                     delta = c('log', 'identity'),
                     theta = c('log'),
                     eta = c('identity'))
    for(par in names(links))
        if(!(links[[par]] %in% par_links[[par]]))
            stop(sprintf('Link %s not allowed for %s', links[[par]], par))
    ## Zwracamy uzupełnioną listę links
    for(par in model_pars[[model]])
        if(!(par %in% names(links)))
            links[[par]] = par_links[[par]][1]
    links
}
## Ok

## parameter dimensionality
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
                 ordinal = c(eta = 1),
                 uvordinal = c(eta = 1, theta = 1))[[model]][par]
        s = rep(s, 2)
    }
    if(any(is.na(s)))
        stop(sprintf('Cannot determine parameter dim: par = %s, model = %s, link = %s, K = %d', par, model, links, K))
    s
}
## Ok
