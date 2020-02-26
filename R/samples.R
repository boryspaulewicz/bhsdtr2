## -*- coding: utf-8 -*-

##' get posterior samples or ml point estimates
##'
##' @export
samples = function(m, par, group = NULL, ...){
    par2 = par.to.linked(par)
    invlink = F
    ## e.g., if par = 'dprim'
    if(par2 != par)
        invlink = T
    if(invlink){
        fun = linkfun(m$links[[par2]])
    }else{
        ## no conversion (e.g., par = 'delta')
        fun = function(x, ...)x
    }
    result = condition.specific.samples(m, par2, group, ...)
    ## par2_size because e.g., gamma_size_ > gamma_size when link = 'parsimonious' or 'twoparameter'
    result_ = array(dim = c(dim(result)[1], m$sdata[[sprintf('%s_size_', par2)]], dim(result)[3]))
    ## i is condition number
    for(i in 1:(dim(result)[3]))
        result_[,,i] = fun(matrix(result[,,i], nrow = dim(result)[1]),
                        m$sdata$K, m$sdata$thresholds_scale)
    dimnames(result_) = dimnames(result)
    dimnames(result_)[2] = list(paste(par, 1:dim(result)[2], sep = '.'))
    class(result_) = c('bhsdtr_samples', class(result_))
    attr(result_, 'data') = attr(result, 'data')
    result_
}

## e.g., dprim -> delta
par.to.linked = function(par)
    if(any(par %in% c('dprim', 'metad', 'sdratio', 'thr', 'mean'))){
        c(dprim = 'delta', metad = 'delta', sdratio = 'theta', thr = 'gamma', mean = 'eta')[par]
    }else{
        par
    }
