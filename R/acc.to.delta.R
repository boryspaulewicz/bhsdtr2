## -*- coding: utf-8 -*-

#' Calculates delta = log(d') given accuracy or accuracy given delta
#'
#' @param acc a vector of average accuracy scores. All the values must
#'     be greater than 0 and lower than 1.
#' @param delta default is NULL, if provided the delta values will be
#'     converted to accuracy
#' @param crit main decision criterion relative to the midpoint
#'     between the evidence distribution means
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
