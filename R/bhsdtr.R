## -*- coding: utf-8 -*-

##' main interface
##' 
##' @export
bhsdtr = function(model_formulae, response_formula, data,
                  links = list(gamma = 'log_distance'), fit_method = 'ml',
                  thresholds_scale = 2, ...){
    fixed = random = list()
    data = as.data.frame(data)
    vs = pars = NULL
    if(length(model_formulae) == 1)
        stop('I need a vector with at least two model formulae!')
    for(i in 1:length(model_formulae)){
        ## fixed effects formula
        ff = lme4::nobars(model_formulae[[i]])
        ## dprim -> delta, sdt ...
        pars = c(pars, as.character(ff[[2]]))
        par = par.to.linked(as.character(ff[[2]]))
        ## ~ fixed effects part ...
        fixed[[par]] = lme4::nobars(ff)[-2]
        ## storing variable names for aggregation
        vs = c(vs, names(get_all_vars(ff[-2], data)))
        ## random effects part, e.g., (1 | id)
        bars = lme4::findbars(model_formulae[[i]])
        if(!is.null(bars)){
            random[[par]] = list()
            for(j in 1:length(bars)){
                ## e.g., ~ duration
                random[[par]][[j]] = c(model.formula = as.formula(paste('~', as.character(bars[[j]])[2])))
                ## e.g., ~ id
                random[[par]][[j]]$group.formula = as.formula(paste('~', as.character(bars[[j]])[3]))
                gmf = model.frame(random[[par]][[j]]$group.formula, data)
                ## e.g., max(id)
                random[[par]][[j]]$group.size = max(fix.index.gaps(gmf[, 1], names(gmf)[1]))
                ## e.g., 'id'
                random[[par]][[j]]$group.name = names(get_all_vars(random[[par]][[j]]$group.formula, data))
                ## just checking
                if(length(random[[par]][[j]]$group.name) > 1)
                    stop(sprintf('More than one grouping variable in %s %s | %s', par, as.character(random[[par]][[j]]$model.formula),
                                 paste(random[[par]][[j]]$group.name, collapse = ' ')))
                ## updating vector of variables for aggregation
                for(k in c('model.formula', 'group.formula'))
                    vs = c(vs, names(get_all_vars(random[[par]][[j]][[k]], data)))
            }
        }
    }
    vs = unique(vs)
    if(setequal(pars, c('dprim', 'thr', 'sdratio'))){
        model = 'uvsdt'
    }else if(setequal(pars, c('dprim', 'thr'))){
        model = 'sdt'
    }else if(setequal(pars, c('metad', 'thr'))){
        model = 'metad'
    }else if(setequal(pars, c('mean', 'thr'))){
        model = 'ordinal'
    }else{
        stop(sprintf('Unknown model. Parameters %s', paste(pars, collapse = ', ')))
    }
    ## Creating resp and stim variables: only resp ~ stim (sdt family)
    ## or resp ~ 1 (ordinal family) allowed
    rsmf = model.frame(response_formula, data)
    rs.vs = names(get_all_vars(response_formula, data))
    if(ncol(rsmf) > 2)
        stop(sprintf('Response formula must be of the form response_variable ~ stimulus_variable or response_variable ~ 1: %s',
                     as.character(response_formula)))
    ## We do not allow for gaps in the response variable e.g., 1, 1,
    ## 2, 4, 5 (no 3s) nor in the stimulus variable: all integers
    ## between 1 and max(stim) or max(resp) have to be present.
    rs = list(resp = fix.index.gaps(rsmf[, 1], names(rs.vs)[1], T))
    if(length(rs.vs) == 2){
        rs$stim = fix.index.gaps(rsmf[, 2], rs.vs[2], T)
    }else{
        ## dummy stimulus variable (response_formula was of the form
        ## resp ~ 1)
        rs$stim = rep(1, nrow(rsmf))
    }
    ## Agreggation
    K = max(rs$resp, na.rm = T)
    links = model.links(model, fixed, random, links, K)
    ## Unique stimulus and response names for aggregation
    rs.names = c(NA, NA)
    ## Only the relevant variables + response and stimulus variables. 
    ## If there is only one variable, we have to make sure data does not
    ## become a vector
    data = data[, vs, drop = F]
    for(i in 1:length(rs.names)){
        ## unique names for resp and stim, _ is safe, some other
        ## choices (e.g., @) result in error
        rs.names[i] = paste(c(names(rs)[i], vs), collapse = '_')
        data[[rs.names[i]]] = rs[[i]]
    }
    res = plyr::ddply(data, unique(c(vs, rs.names[names(rs) == 'stim'])),
                      ## We are adding the 1:K vector to make sure
                      ## that the counts for every k are here, the -1
                      ## term corrects for this
                      function(df)table(c(df[[rs.names[names(rs) == 'resp']]], 1:K)) - 1)
    counts = res[, c((ncol(res) - K + 1):ncol(res))]
    ## aggregated data object
    adata = list(data = res[, setdiff(vs, rs.names), drop = F])
    adata$stimulus = res[[rs.names[names(rs) == 'stim']]]
    ## Creatinng data structures and model code for stan
    sdata = list(N = nrow(adata$data), K = K, Kb2 = round(K / 2), PRINT = 0,
                 unbiased = fix.stan.dim(unbiased(K)), thresholds_scale = thresholds_scale,
                 counts = counts)
    ## in SDT models stim_sign = -1, 1 (just an ugly hack)
    sdata$stim_sign = 2 * as.numeric(as.factor(as.character(adata$stimulus))) - 3
    ## Fixed effects' model matrices
    for(par in names(fixed)){
        v = sprintf('X_%s', par)
        sdata[[v]] = model.matrix(fixed[[par]], adata$data)
        sdata[[sprintf('X_%s_ncol', par)]] = ncol(sdata[[v]])
        sdata[[sprintf('%s_is_fixed', par)]] = sdata[[sprintf('%s_fixed_value', par)]] =
            matrix(0, nrow = par.size(par, model, links, K)[1], ncol = ncol(sdata[[v]]))
        ## priors
        for(prior.par in c('fixed_mu', 'fixed_sd'))
            sdata[[sprintf('%s_prior_%s', par, prior.par)]] = default_prior(par, ncol(sdata[[v]]), prior.par, model, links, K)
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
                g[[i]] = fix.index.gaps(mf[, 1], names(mf)[1])
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
            sdata[[sprintf('%s_prior_scale_%d', par, i)]] = default_prior(par, Z_ncol[i], 'random_scale', model, links, K)
        }
    }
    ## Matrix dimensions
    for(par in unique(c(names(fixed), names(random)))){
        sdata[[sprintf('%s_size', par)]] = par.size(par, model, links, K)[1]
        sdata[[sprintf('%s_size_', par)]] = par.size(par, model, links, K)[2]
    }
    if(is.null(sdata$delta_size)){
        sdata$dprim_size = 1
    }else{
        sdata$dprim_size = sdata$delta_size
    }
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

model.links = function(model, fixed, random, links, K){
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
