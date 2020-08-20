##' Fit a cumulative ordinal model with varying ordered thresholds and fixed latent mean and sd
##'
##' @param formula a model formula of the form resp ~ f1 * f2 + (f1 |
##'     id), where resp is the name of the ordinal response
##'     variable. This will be transformed to a vector of two model
##'     formulae and a response - stimulus formula for the bhsdtr
##'     function: bhsdtr(c(mean ~ 1, thr ~ f1 * f2 + (f1 | id)), resp
##'     ~ 1, ...)
##' @param ... additional arguments passed to the bhsdtr function. See
##'     the bhsdtr documentation for details.
##' @return a bhsdtr model object
##' @export
cumulative = function(formula, ...){
    resp = as.character(formula[[2]])
    theta.form = as.formula(paste(c('thr', as.character(formula[-2])), collapse = ' '))
    bhsdtr(c(mean ~ 1, theta.form), as.formula(sprintf('%s ~ 1', resp)), ...)
}
