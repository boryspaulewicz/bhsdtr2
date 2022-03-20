## -*- coding: utf-8 -*-

##' get posterior samples or jmap point estimates per condition
##'
##' @export
condition.specific.samples = function(m, par, group = NULL, method = NULL, include = NULL){
    if(is.null(method)){
        if(!is.null(m$stanfit)){
            method = 'stan'
        }else if(!is.null(m$jmapfit)){
            method = 'jmap'
        }else{
            stop('This model was not fitted')
        }
    }
    if(method == 'stan'){
        slot = 'stanfit'
    }else if(method == 'jmap'){
        slot = 'jmapfit'
    }else if(method == 'prior'){
        slot = 'stanfit.prior'
        if(is.null(m$stanfit.prior))
            stop('Prior was not sampled for this model')
    }
    fixed.formula = m$fixed[[par]]
    vnames = unique(c(names(get_all_vars(fixed.formula, m$adata$data)), include))
    if(!is.null(group)){
        group.name = m$random[[par]][[group]][['group.name']]
        group.index = m$sdata[[sprintf('%s_group_%d', par, group)]]
    }else{
        group.name = NULL
    }
    vnames = sort(unique(c(vnames, group.name)))
    data = m$adata$data[, vnames, drop = F]
    if(!is.null(group)){
        group.to.index = unique(data.frame(original = data[[group.name]], index = group.index))
        rownames(group.to.index) = as.character(group.to.index$original)
        ## data[[group.name]] = group.index
    }
    data = unique(data[, vnames, drop = F])
    ## if e.g., dprim ~ 1
    if(nrow(data) == 0)
        data = data.frame(x = '')
    ## colnames describe unique combinations of predictors
    condition.names = apply(data, 1, function(x)paste(x, collapse = ':'))
    X = model.matrix(fixed.formula, data)
    if(('stan' %in% method) & !is.null(m[[slot]])){
        ## samples.fixef = extract(m[[slot]])[[sprintf('%s_fixed', par)]]
        samples.fixef = merged.extract(m, par)
    }else if(('jmap' %in% method) & !is.null(m[[slot]])){
        samples.fixef = array(m[[slot]]$par[grep(sprintf('%s_fixed\\[', par), names(m[[slot]]$par))],
                   dim = c(1, m$sdata[[sprintf('%s_size', par)]], ncol(X)))
    }else{
        stop(sprintf('This model was not fitted using method %s', paste(method, collapse = ' nor ')))
    }
    if(!is.null(group)){
        random.formula = m$random[[par]][[group]][['model.formula']]
        Z = model.matrix(random.formula, data)
        if(('stan' %in% method) & !is.null(m[[slot]])){
            ## samples.ranef = extract(m[[slot]])[[sprintf('%s_random_%d', par, group)]]
            samples.ranef = merged.extract(m, par, group) 
        }else if('jmap' %in% method){
            samples.ranef = array(m[[slot]]$par[grep(sprintf('%s_random_%d\\[', par, group), names(m[[slot]]$par))],
                       dim = c(1, m$random[[par]][[group]]$group.size, m$sdata[[sprintf('%s_size', par)]],
                               ncol(X)))
        }
    }
    ## number of samples, par.size, number of conditions, where condition may include group level
    result = result.ranef =
        array(dim = c(dim(samples.fixef)[1], dim(samples.fixef)[2], nrow(data)))
    for(s in 1:(dim(samples.fixef)[1]))
        result[s,,] = samples.fixef[s,,] %*% t(X)
    if(!is.null(group)){
        group.index = group.to.index[as.character(data[[group.name]]), 'index']
        for(s in 1:(dim(samples.fixef)[1]))
            for(con in 1:nrow(data)){
                ##! an ugly hack
                if(dim(samples.ranef)[4] > 1){
                    result.ranef[s,,con] =
                        samples.ranef[s, group.index[con],,] %*% t(Z[con,,drop = F])
                }else{
                    result.ranef[s,,con] =
                        ## samples.ranef[s, data[[group.name]][con],,,drop = F] %*% t(Z[con,,drop = F])
                        samples.ranef[s, group.index[con],,,drop = F] %*% t(Z[con,,drop = F])
                }
            }
        if(m$links[[par]] == 'id_log'){
            result = result * exp(result.ranef)
        }else{
            result = result + result.ranef
        }
    }
    dimnames(result)[3] = list(condition.names)
    dimnames(result)[2] = list(paste(par, 1:dim(result)[2], sep = '.'))
    attr(result, 'data') = data
    class(result) = c('bhsdtr_samples', class(result))
    result
}
