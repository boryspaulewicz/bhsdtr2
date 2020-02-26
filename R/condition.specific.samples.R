## -*- coding: utf-8 -*-

##' get posterior samples or ml point estimates per condition
##'
##' @export
condition.specific.samples = function(m, par, group = NULL, fit = c('stan', 'ml'), include.vars = NULL){
    fixed.formula = m$fixed[[par]]
    vnames = unique(c(names(get_all_vars(fixed.formula, m$adata$data)), include.vars))
    if(!is.null(group)){
        group.name = m$random[[par]][[group]][['group.name']]
        group.index = m$sdata[[sprintf('%s_group_%d', par, group)]]
    }else{
        group.name = NULL
    }
    data = m$adata$data[, unique(c(vnames, group.name)), drop = F]
    if(!is.null(group))
        data[[group.name]] = group.index
    data = unique(data[, unique(c(group.name, vnames)), drop = F])
    ## if par ~ 1
    if(nrow(data) == 0)
        data = data.frame(x = '')
    ## colnames describe unique conditions
    condition.names = apply(data, 1, function(x)paste(x, collapse = ':'))
    X = model.matrix(fixed.formula, data)
    if(('stan' %in% fit) & !is.null(m$stanfit)){
        samples.fixef = extract(m$stanfit)[[sprintf('%s_fixed', par)]]
    }else if(('ml' %in% fit) & !is.null(m$mlfit)){
        samples.fixef = array(m$mlfit$par[grep(sprintf('%s_fixed\\[', par), names(m$mlfit$par))],
                   dim = c(1, m$sdata[[sprintf('%s_size', par)]], ncol(X)))
    }else{
        stop(sprintf('This model was not fitted using method %s', paste(fit, collapse = ' nor ')))
    }
    if(!is.null(group)){
        random.formula = m$random[[par]][[group]][['model.formula']]
        Z = model.matrix(random.formula, data)
        if(('stan' %in% fit) & !is.null(m$stanfit)){
            samples.ranef = extract(m$stanfit)[[sprintf('%s_random_%d', par, group)]]
        }else if('ml' %in% fit){
            samples.ranef = array(m$mlfit$par[grep(sprintf('%s_random_%d\\[', par, group), names(m$mlfit$par))],
                       dim = c(1, m$random[[par]][[group]]$group.size, m$sdata[[sprintf('%s_size', par)]],
                               ncol(X)))
        }
    }
    ## number of samples, par.size, number of conditions, where condition may include group level
    result = array(dim = c(dim(samples.fixef)[1], dim(samples.fixef)[2], nrow(data)))
    ## result = array(dim = c(dim(samples.fixef)[1:2], nrow(data)))
    for(s in 1:(dim(samples.fixef)[1])){
        result[s,,] = samples.fixef[s,,] %*% t(X)
        if(!is.null(group)){
            ## adding the random effects
            for(con in 1:nrow(data)){
                result[s,,con] = result[s,,con] +
                    samples.ranef[s, data[[group.name]][con],,] %*% t(Z[con, ])
            }
        }
    }
    dimnames(result)[3] = list(condition.names)
    dimnames(result)[2] = list(paste(par, 1:dim(result)[2], sep = '.'))
    attr(result, 'data') = data
    class(result) = c('bhsdtr_samples', class(result))
    result
}
