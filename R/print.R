## -*- coding: utf-8 -*-

##' print a bhsdtr model object
##'
##' @export
print.bhsdtr_model = function(x, ...){
    cat('bhsdtr model object\n')
    cat(sprintf('model type: %s\n', x$model))
    cat(sprintf('K: %d\n', x$K))
    cat('link functions:\n')
    for(par in sort(names(x$links)))
        cat(sprintf(' %s: %s\n', par, x$links[[par]]))
    for(v in sort(names(x$sdata)[grep('[a-z]+_size', names(x$sdata))]))
        cat(sprintf('%s: %d\n', v, x$sdata[[v]]))
    cat('fixed effects:\n')
    for(par in sort(names(x$fixed)))
        cat(sprintf(' %s %s\n', par, paste(as.character(x$fixed[[par]]), collapse = ' ')))
    cat('random effects:\n')
    for(par in sort(names(x$random)))
        for(i in 1:length(x$random[[par]]))
            cat(sprintf(' %s %s | %s\n', par, paste(as.character(x$random[[par]][[i]]$model.formula), collapse = ' '),
                        x$random[[par]][[i]]$group.name))
    cat(sprintf('variables: %s\n', paste(sort(names(x$adata$data)), collapse = ' ')))
    cat(sprintf('original data size: %d\n', x$data_size))
    cat(sprintf('aggregated data size: %d\n', nrow(x$adata$data)))
    if(!is.null(x$mlfit))
        cat('mlfit\n')
    if(!is.null(x$stanfit)){
        cat('preparing stan summary statistics...\n')
        print(x$stanfit, probs = c(.025, .957),
              pars = x$pars[-c(grep('counts_new', x$pars), grep('_random_', x$pars), grep('corr_', x$pars))])
    }
    invisible(x)
}
## Ok
