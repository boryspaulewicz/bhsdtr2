## -*- coding: utf-8 -*-

##' Link function by name
##'
##' @export
linkfun = function(link){
    res = NULL
    if(link == 'log')
        res = function(x, ...)exp(x)
    if(link == 'identity')
        res = function(x, ...)x
    if(link == 'log_ratio')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            criteria[, K / 2] = samples[, K / 2];
            if(K > 2){
                ## spread
                criteria[, K / 2 + 1] = criteria[, K / 2] + exp(samples[, K / 2 + 1]);
                ## symmetry
                criteria[, K / 2 - 1] = criteria[, K / 2] - exp(samples[, K / 2 - 1]) * (criteria[, K / 2 + 1] - criteria[, K / 2]);
                if(K > 4){
                    for(k in 1:(K / 2 - 2)){
                        ## upper consistency
                        criteria[, K / 2 + k + 1] = criteria[, K / 2 + k] + exp(samples[, K / 2 + k + 1]) * (criteria[, K / 2 + 1] - criteria[, K / 2]);
                        ## lower consistency
                        criteria[, K / 2 - k - 1] = criteria[, K / 2 - k] - exp(samples[, K / 2 - k - 1]) * (criteria[, K / 2] - criteria[, K / 2 - 1]);
                    }
                }
            }
            criteria
        }
    if(link == 'log_distance')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            criteria[, K / 2] = samples[, K / 2];
            if(K > 2){
                for(k in 1:(K / 2 - 1)){
                    criteria[, K / 2 + k] = criteria[, K / 2 + k - 1] + exp(samples[, K / 2 + k]);
                    criteria[, K / 2 - k] = criteria[, K / 2 - k + 1] - exp(samples[, K / 2 - k]);
                }
            }
            criteria
        }
    if(link == 'parsimonious')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            unb = unbiased(K)
            for(k in 1:(K - 1))
                criteria[, k] = samples[, 1] + exp(samples[, 2]) * unb[k]
            criteria
        }
    if(link == 'twoparameter')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            criteria[, K / 2] = samples[, 1];
            if(K > 2){
                for(k in 1:(K / 2 - 1)){
                    criteria[, K / 2 + k] = criteria[, K / 2 + k - 1] + exp(samples[, 2]);
                    criteria[, K / 2 - k] = criteria[, K / 2 - k + 1] - exp(samples[, 2]);
                }
            }
            criteria
        }
    if(link == 'softmax')
        res = function(samples, K, thresholds_scale)
            t(apply(exp(cbind(samples, 0)), 1,
                    function(x) thresholds_scale * stats::qnorm(cumsum(x/sum(x))[-length(x)])))
    res
}
