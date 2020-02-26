## -*- coding: utf-8 -*-

##' get posterior samples or ml point estimates
##'
##' @export
samples = function(m, par, group = NULL){
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
    res = condition.specific.samples(m, par2, group)
    if(is.null(group)){
        ## e.g., gamma_size_ > gamma_size when link = 'parsimonious' or 'twoparameter'
        res_ = array(dim = c(dim(res)[1], m$sdata[[sprintf('%s_size_', par2)]], dim(res)[3]))
        ## i is condition number
        for(i in 1:(dim(res)[3]))
            res_[,,i] = fun(matrix(res[,,i], nrow = dim(res)[1]),
                           m$sdata$K, m$sdata$thresholds_scale)
    }else{
        res_ = array(dim = c(dim(res)[1:2], m$sdata[[sprintf('%s_size_', par2)]], dim(res)[4]))
        for(g in 1:(dim(res)[2]))
            for(i in 1:(dim(res)[4]))
                res_[,g,,i] = fun(matrix(res[,g,,i], nrow = dim(res)[1]),
                                  m$sdata$K, m$sdata$thresholds_scale)
    }
    dimnames(res_) = dimnames(res)
    res_
}

## e.g., dprim -> delta
par.to.linked = function(par)
    if(any(par %in% c('dprim', 'metad', 'sdratio', 'thr', 'mean'))){
        c(dprim = 'delta', metad = 'delta', sdratio = 'theta', thr = 'gamma', mean = 'eta')[par]
    }else{
        par
    }
