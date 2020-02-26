## -*- coding: utf-8 -*-

#' Calculates delta = log(d') given accuracy or accuracy given delta
#' assuming no bias and p(stim = 1) = 0.5.
#' @param acc a vector of average accuracy scores. All the values must
#'     be greater than 0 and lower than 1.
#' @param delta default is NULL, if provided the delta values will be
#'     converted to accuracy of an unbiased observer
#' @return delta = log(d') or expected accuracy
#' @export
acc.to.delta = function(acc = .75, delta = NULL, crit = 0){
    if(any((acc <= 0) || (acc >= 1)))
        stop('Some accuracy scores were to extreme')
    if(is.null(delta)){
        log(2 * stats::qnorm(acc))
    }else{
        (pnorm(exp(delta) / 2 - crit) + pnorm(exp(delta) / 2 + crit)) / 2
    }
}
