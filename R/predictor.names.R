##' Extract the names of the predictor variables
##'
##' @param model a bhsdtr model object
##' @param ... e.g., 'delta', 'gamma'
##' @return vector of names of variables in the original data frame
##' @export
predictor.names = function(model, ...){
    vnames = NULL
    for(v in list(...))
        vnames = c(vnames, names(get_all_vars(model$fixed[[v]], model$data)))
    vnames
}
