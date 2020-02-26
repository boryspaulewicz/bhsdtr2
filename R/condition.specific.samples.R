## -*- coding: utf-8 -*-

##' get posterior samples or ml point estimates per condition
##'
##' @export
condition.specific.samples = function(m, par, group = NULL){
    ff = m$fixed[[par]]
    df = unique(get_all_vars(ff, m$adata$data))
    ## if par ~ 1
    if(nrow(df) == 0)
        df = data.frame(x = 'x')
    ## colnames describe conditions
    cnames = apply(df, 1, function(x)paste(x, collapse = ':'))
    ## Macierz efektów ustalonych. Ogólnie, macierz będzie kwadratowa
    ## tylko dla samych czynników
    X = model.matrix(ff, df)
    if(!is.null(m$stanfit)){
        sf = extract(m$stanfit)[[sprintf('%s_fixed', par)]]
    }else if(!is.null(m$mlfit)){
        sf = array(m$mlfit$par[grep(sprintf('%s_fixed\\[', par), names(m$mlfit$par))],
                   dim = c(1, m$sdata[[sprintf('%s_size', par)]], ncol(X)))
    }else{
        stop('This model was not fitted')
    }
    if(!is.null(group)){
        fr = m$random[[par]][[group]][['model.formula']]
        Z = model.matrix(fr, df)
        if(!is.null(m$stanfit)){
            sr = extract(m$stanfit)[[sprintf('%s_random_%d', par, group)]]
        }else if(!is.null(m$mlfit)){
            sr = array(m$mlfit$par[grep(sprintf('%s_random_%d\\[', par, group), names(m$mlfit$par))],
                       dim = c(1, m$random[[par]][[group]]$group.size, m$sdata[[sprintf('%s_size', par)]],
                               ncol(X)))
        }else{
            stop('This model was not fitted')
        }
        ## number of samples, number of groups, par.size, number of conditions
        res = array(dim = c(dim(sf)[1], dim(sr)[2], dim(sf)[2], nrow(df)))
    }else{
        ## number of samples, par.size, number of conditions
        res = array(dim = c(dim(sf)[1], dim(sf)[2], nrow(df)))
    }
    ## res = array(dim = c(dim(sf)[1:2], nrow(df)))
    for(s in 1:(dim(sf)[1])){
        if(is.null(group)){
            res[s,,] = sf[s,,] %*% t(X)
        }else{
            for(g in 1:m$random$gamma[[group]]$group.size){
                res[s,g,,] = sf[s,,] %*% t(X)
                res[s,g,,] = sr[s,g,,] %*% t(Z)
            }
        }
    }
    if(is.null(group)){
        dimnames(res)[3] = list(cnames)
    }else{
        dimnames(res)[4] = list(cnames)
    }
    res
}
