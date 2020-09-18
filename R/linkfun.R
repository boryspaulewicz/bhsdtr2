## -*- coding: utf-8 -*-

##' Link function by name
##'
##' @export
linkfun = function(link){res = NULL
    if(link == 'identity') res = function(samples, ...)samples
    ## if id_log, then condition.specific.samples takes care of the log part for random effects
    if(link == 'id_log') res = function(samples, ...)samples
    if(link == 'log') res = function(samples, ...)exp(samples)
    if(link == 'logit') res = function(samples, ...)binomial()$linkinv((samples))
    if(link == 'log_logit')
        res = function(samples, ...)
            cbind(exp(samples[, 1]), binomial()$linkinv((samples[, 2])))
    if(link == 'log_ratio')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            Kb2 = round(K / 2)
            criteria[, Kb2] = samples[, Kb2];
            if(K > 2){
                ## spread
                criteria[, Kb2 + 1] = criteria[, Kb2] + exp(samples[, Kb2 + 1]);
                ## symmetry
                criteria[, Kb2 - 1] = criteria[, Kb2] - exp(samples[, Kb2 - 1]) * (criteria[, Kb2 + 1] - criteria[, Kb2]);
                if(K > 4){
                    for(k in 1:(Kb2 - 2)){
                        ## upper consistency
                        criteria[, Kb2 + k + 1] = criteria[, Kb2 + k] + exp(samples[, Kb2 + k + 1]) * (criteria[, Kb2 + 1] - criteria[, Kb2]);
                        ## lower consistency
                        criteria[, Kb2 - k - 1] = criteria[, Kb2 - k] - exp(samples[, Kb2 - k - 1]) * (criteria[, Kb2] - criteria[, Kb2 - 1]);
                    }
                }
            }
            criteria
        }
    if(link == 'log_distance')
        res = function(samples, K, ...){
            criteria = matrix(nrow = nrow(samples), ncol = K - 1)
            Kb2 = round(K / 2)
            criteria[, Kb2] = samples[, Kb2];
            if(K > 2){
                for(k in 1:(Kb2 - 1)){
                    criteria[, Kb2 + k] = criteria[, Kb2 + k - 1] + exp(samples[, Kb2 + k]);
                    criteria[, Kb2 - k] = criteria[, Kb2 - k + 1] - exp(samples[, Kb2 - k]);
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
            Kb2 = round(K / 2)
            criteria[, round(K / 2)] = samples[, 1];
            if(K > 2){
                for(k in 1:(Kb2 - 1)){
                    criteria[, Kb2 + k] = criteria[, Kb2 + k - 1] + exp(samples[, 2]);
                    criteria[, Kb2 - k] = criteria[, Kb2 - k + 1] - exp(samples[, 2]);
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
