## -*- coding: utf-8 -*-

##' print a bhsdtr model object
##'
##' @export
print.bhsdtr_model = function(x, ...){
    cat(sprintf('bhsdtr model of type %s, number of response categories K: %d\n', x$model, x$sdata$K))
    cat('link functions:\n')
    str = NULL
    for(par in sort(names(x$links)))
        str = c(str, sprintf(' %s: %s', par, x$links[[par]]))
    cat(paste(str, collapse = ', '))
    cat('\n')
    ## for(v in sort(names(x$sdata)[grep('[a-z]+_size', names(x$sdata))]))
    ##     cat(sprintf('%s: %d\n', v, x$sdata[[v]]))
    cat('fixed effects:\n')
    str = NULL
    for(par in sort(names(x$fixed)))
        str = c(str, sprintf(' %s %s', par, paste(as.character(x$fixed[[par]]), collapse = ' ')))
    cat(paste(str, collapse = ', '))
    cat('\n')
    cat('random effects:\n')
    str = NULL
    for(par in sort(names(x$random)))
        for(i in 1:length(x$random[[par]]))
            str = c(str, sprintf(' %s %s | %s', par, paste(as.character(x$random[[par]][[i]]$model.formula), collapse = ' '),
                        x$random[[par]][[i]]$group.name))
    cat(paste(str, collapse = ', '))
    cat('\n')
    cat(sprintf('variables: %s, ', paste(sort(names(x$adata$data)), collapse = ' ')))
    cat(sprintf('data size original: %d, aggregated: %d\n', nrow(x$data), nrow(x$adata$data)))
    str = NULL
    if(!is.null(x$jmapfit))
        str = c(str, 'jmap')
    if(!is.null(x$stanfit))
        str = c(str, 'stan')
    if(!is.null(str))
        cat(sprintf('fitted using method: %s\n', paste(str, collapse = ', ')))
    invisible(x)
}

##' print bhsdtr condition-specific samples or jmap point estimates
##'
##' @export
print.bhsdtr_samples = function(x, digits = 2, probs = c(.025, .975), ...){
    cat(sprintf('samples: %d, estimates rounded to %d decimal places\n', dim(x)[1], digits))
        tbl = matrix(x[1,,], nrow = dim(x)[2], ncol = dim(x)[3])
        for(d in 1:(dim(x)[2]))
            tbl[d,] = apply(x[,d,, drop = F], 3, mean)
        tbl = as.data.frame(formatC(round(tbl, digits), format = 'f', digits = digits))
        rownames(tbl) = dimnames(x)[[2]]
        names(tbl) = dimnames(x)[[3]]
        print(as.data.frame(t(tbl)))
    invisible(x)
}
