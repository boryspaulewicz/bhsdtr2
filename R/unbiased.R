## -*- coding: utf-8 -*-

##' "unbiased" thresholds for the parsimonious link function
##' 
##' @export
unbiased = function(K){
    log((1:K / K) / (1 - 1:K / K))[1:(K - 1)]
}
