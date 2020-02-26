## -*- coding: utf-8 -*-

##' summary of condition specific samples
##'
##' @export
summary.bhsdtr_samples = function(x, digits = 1, probs = c(.025, .975), ...){
        tbl = as.data.frame(matrix(x[1,,], nrow = dim(x)[2], ncol = dim(x)[3]))
        format.fun = function(x)paste(formatC(round(x, digits), format = 'f', digits = digits,  width = digits + 3),
                                      collapse = ' ')
        if(dim(x)[1] == 1){
            sumfun = function(x)format.fun(mean(x))
        }else{
            sumfun = function(x)format.fun(c(mean(x), quantile(x, probs))[c(2, 1, 3)])
        }
        for(d in 1:(dim(x)[2]))
            tbl[d,] = apply(x[,d,, drop = F], 3, sumfun)
        rownames(tbl) = dimnames(x)[[2]]
        names(tbl) = dimnames(x)[[3]]
        tbl
}
