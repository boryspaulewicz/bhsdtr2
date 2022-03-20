## -*- coding: utf-8 -*-

##' plot a fitted bhsdtr model
##'
##' @export
plot.bhsdtr_model = function(x, vs = NULL, type = 'response', alpha = .05, bw = FALSE, verbose = T, fit = 'default', ...){
    ## only ~ 1, so no variables in the aggregated data object
    if(!is.null(vs))
        if(all(vs == 1)){
            x$adata$data[paste(names(x$adata$data), collapse = ':')] = 1
            vs = tail(names(x$adata$data), 1)
        }
    if(ncol(x$adata$data) == 0)
        x$adata$data = data.frame(x = rep(' ', nrow(x$adata$data)))
    if(fit == 'default'){
        if(!is.null(x$stanfit)){
            fit = 'stan'
        }else if(!is.null(x$jmapfit)){
            fit = 'jmap'
        }else{
            stop('This model was not fitted')
        }
    }
    if(fit == 'jmap'){
        p = pointest.plot(x$jmapfit, x$code, x$model, x$adata, x$sdata, vs, bw)
    }else{
        if(is.null(vs)){
            vs = names(x$adata$data)
            for(rpar in x$random)
                for(rg in rpar)
                    vs = setdiff(vs, rg$group.name)
        }
        s = as.data.frame(x$stanfit)
        cnt_new = t(s[,grep('counts_new', names(s))])
        rm(s)
        ## Tworzymy zbiór danych rozwinięty pionowo
        df = cbind(x$adata$data, x$adata$stimulus)
        names(df)[ncol(df)] = 'stimulus'
        stim = ncol(df)
        K = ncol(x$sdata$counts)
        df = cbind(plyr::ddply(df, names(df), function(x)data.frame(1:K)), as.vector(t(x$sdata$counts)))
        df = cbind(df, rep(1:nrow(x$adata$data), each = ncol(x$sdata$counts)))
        ## Indeksy ważnych kolumn
        resp = ncol(df) - 2
        cnt = ncol(df) - 1
        obs = ncol(df)
        ## Kolejność taka, jak w próbkach ze Stan-a
        df = df[order(df[[resp]], df[[obs]]),]
        ## Agregacja obserwacji i próbek dla dowolnego wykresu
        if(length(vs) > 0){
            f = as.factor(df[[vs[1]]])
            if(length(vs) > 1)
                for(v in vs[-1])
                    f = f:as.factor(df[[v]])
        }else{
            f = as.factor(rep('', nrow(df)))
        }
        dfa = stats::aggregate(df[[cnt]] ~ df[[resp]] + df[[stim]] + f, FUN = sum)
        cnt_new_a = matrix(nrow = ncol(cnt_new), ncol = nrow(dfa))
        if(verbose)pb = utils::txtProgressBar(min = 1, max = ncol(cnt_new), style = 3)
        if(verbose)print('Aggregating posterior samples...')
        for(r in 1:ncol(cnt_new)){
            cnt_new_a[r,] = stats::aggregate(cnt_new[,r] ~ df[[resp]] + df[[stim]] + f, FUN = sum)[[4]]
            if(verbose)utils::setTxtProgressBar(pb, r)
        }
        if(verbose)close(pb)
        rm(cnt_new)
        names(dfa) = c('response', 'stimulus', 'f', 'count')
        dfa$stimulus = as.factor(dfa$stimulus)
        dfa$n = plyr::ddply(dfa, c('f', 'stimulus'), function(x)data.frame(n = rep(sum(x$count), nrow(x))))[[3]]
        dfa$i = 1:nrow(dfa)
        if((type == 'roc') & (x$model != 'ordinal')){
            ## Wykres ROC: wyliczamy p(Hit) ~ p(FA)
            if(verbose)print('Calculating ROC curves...')
            dfroc = plyr::ddply(dfa, c('stimulus', 'f'),
                                function(x){
                                    cumfr_new = apply(t(cnt_new_a[,x$i]), 2, function(v)rev(cumsum(rev(v / x$n))))
                                    data.frame(response = x$response,
                                               cumfr = rev(cumsum(rev(x$count / x$n))),
                                               cumfr.fit = apply(cumfr_new, 1, mean),
                                               cumfr.lo = apply(cumfr_new, 1, function(x)stats::quantile(x, alpha / 2)),
                                               cumfr.hi = apply(cumfr_new, 1, function(x)stats::quantile(x, 1 - alpha / 2)))
                                }) ## .progress = 'text'
            rm(dfa)
            ## dfroc$in.pi = as.numeric((dfroc$cumfr >= dfroc$cumfr.lo) && (dfroc$cumfr <= dfroc$cumfr.hi))
            ## dfroc$in.pi[dfroc$in.pi == 0] = .5
            dfrocs = dfroc[dfroc$stimulus == '1',]
            dfrocs$stim2 = dfroc[dfroc$stimulus == '2',]
            rm(dfroc)
            p = ggplot2::ggplot(dfrocs, ggplot2::aes(cumfr, stim2$cumfr)) +
                ggplot2::geom_line(ggplot2::aes(x = cumfr.fit, y = stim2$cumfr.fit), lty = 2) +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = stim2$cumfr.lo, ymax = stim2$cumfr.hi, x = cumfr.fit), width = 0.02) +
                ggplot2::geom_errorbarh(ggplot2::aes(xmin = cumfr.lo, xmax = cumfr.hi, y = stim2$cumfr.fit), height = 0.02) +
                ggplot2::geom_point() +
                ggplot2::labs(x = 'p(F)', y = 'p(H)') +
                ggplot2::coord_fixed() +
                ggplot2::facet_wrap(~f)
            if(bw){
                p + theme_minimalist()
            }else{
                p
            }
        }else{
            ## Wykres rozkładów odpowiedzi
            dfa[,c('count.lo', 'count.hi', 'count.fit')] = cbind(apply(cnt_new_a, 2, function(x)stats::quantile(x, alpha / 2)),
                                                                 apply(cnt_new_a, 2, function(x)stats::quantile(x, 1 - alpha / 2)),
                                                                 apply(cnt_new_a, 2, mean))
            dfa$response = dfa$response + (as.numeric(dfa$stimulus) - 1.5) / 4
            if(bw){
                p = ggplot2::ggplot(dfa, ggplot2::aes(response, count / n), group = stimulus) +
                    ggplot2::geom_errorbar(ggplot2::aes(ymin = count.lo / n, ymax = count.hi / n, lty = stimulus), width = 0.2) +
                    ggplot2::geom_line(ggplot2::aes(y = count.fit / n, lty = stimulus)) +
                    ggplot2::geom_point(ggplot2::aes(pch = stimulus)) +
                    ggplot2::labs(x = 'Response', y = 'Frequency', pch = 'Stimulus', lty = 'Stimulus') +
                    ggplot2::facet_wrap(~f) + theme_minimalist()
            }else{
                p = ggplot2::ggplot(dfa, ggplot2::aes(response, count / n, color = stimulus), group = stimulus) +
                    ggplot2::geom_errorbar(ggplot2::aes(ymin = count.lo / n, ymax = count.hi / n), width = 0.2) +
                    ggplot2::geom_line(ggplot2::aes(y = count.fit / n)) +
                    ggplot2::geom_point(ggplot2::aes(pch = stimulus)) +
                    ggplot2::labs(x = 'Response', y = 'Frequency', color = 'Stimulus', pch = 'Stimulus') +
                    ggplot2::facet_wrap(~f)
            }
            p
        }
    }
    p
}

pointest.plot = function(jmapfit, model_code, model, adata, sdata, vs = NULL, bw = FALSE){
    pointest = matrix(jmapfit$par[grep('multinomial_p', names(jmapfit$par))], nrow = max(1, nrow(adata$data)))
    obs = t(apply(sdata$counts, 1, function(x)(x / sum(x))))
    sdt_models = c('sdt', 'uvsdt', 'metad', 'dpsdtcor', 'dpsdt')
    if(model %in% sdt_models)
        stim = as.factor(adata$stimulus)
    if(is.null(vs)){
        f = as.factor(adata$data[,1])
        if(ncol(adata$data) > 1)
            for(i in 2:ncol(adata$data))
                f = f:as.factor(adata$data[,i])
    }else{
        f = as.factor(adata$data[[vs[1]]])
        if(length(vs) > 1)
            for(v in 2:length(vs))
                f = f:as.factor(adata$data[[vs[v]]])
    }
    df = data.frame(f = rep(f, ncol(pointest)), th = rep(1:ncol(pointest), each = length(f)),
                    pointest = as.vector(pointest), obs = as.vector(obs))
    aggregation.vs = c('th', 'f')
    if(model %in% sdt_models){
        df$stim = rep(stim, ncol(pointest))
        aggregation.vs = c(aggregation.vs, 'stim')
    }
    adf = aggregate(df$obs, df[, aggregation.vs], mean)
    adf$fit = aggregate(df$pointest, df[, aggregation.vs], mean)$x
    title = sprintf('Lines = model: %s, points = observed response distribution', model)
    yl = 'p(Response)'
    if(model %in% sdt_models){
        xl = 'Response = Decision + Rating'
        p = ggplot2::ggplot(adf, ggplot2::aes(th, x, group = stim, color = stim)) +
            ggplot2::geom_line(ggplot2::aes(y = fit)) +
            ggplot2::geom_point() +
            ggplot2::labs(title = title, color = 'Stimulus') +
            ggplot2::xlab(xl) +
            ggplot2::ylab(yl) +
            ggplot2::facet_wrap(~ f)
    }else{
        xl = 'Response = Rating'
        p = ggplot2::ggplot(adf, ggplot2::aes(th, x)) +
            ggplot2::geom_line(ggplot2::aes(y = fit)) +
            ggplot2::geom_point() +
            ggplot2::labs(title = title) +
            ggplot2::xlab(xl) +
            ggplot2::ylab(yl) +
            ggplot2::facet_wrap(~ f)
    }
    if(bw){
        p + theme_minimalist()
    }else{
        p
    }
}

