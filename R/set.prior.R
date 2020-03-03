##' Set the parameters of the prior distributions
##' 
##' @param model a bhsdtr model object
##' @param ... e.g., delta_prior_fixed_mu = log(1.5),
##'     delta_prior_fixed_sd = 2. Note that the prior parameters are
##'     internally represented by matrices of dimension DxC, where D
##'     is the dimensionality of the parameter (delta has 2 in the
##'     meta-d' model, otherwise it has 1, ) and C is the number of
##'     columns in the fixed (<par>_prior_fixed_mu,
##'     <par>_prior_fixed_sd) or the random (<par>prior_scale_<g>)
##'     effects' model matrix. You can provide scalars, vectors, or
##'     matrices. A vector will fill the prior parameter matrix in
##'     column major order.
##' @return a bhsdtr model object with updated priors
##' @export
set.prior = function(model, ...){
    args = list(...)
    for(arg in names(args)){
        if(length(grep('prior', arg)) == 0)
            stop(sprintf('%s is not a prior parameter', arg))
        if(is.null(model$sdata[[arg]])){
            stop(sprintf('Prior parameter %s does not exist in this model', arg))
        }else{
            model$sdata[[arg]] = matrix(args[[arg]], nrow = nrow(model$sdata[[arg]]), ncol = ncol(model$sdata[[arg]]))
        }
    }
    model
}
