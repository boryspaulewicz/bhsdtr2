## -*- coding: utf-8 -*-

##' plot a fitted bhsdtr model
##'
##' @export
plot.bhsdtr_model = function(x, vs = NULL, type = 'response', alpha = .05, bw = FALSE, verbose = T, ...){
    ## only ~ 1, so no variables in the aggregated data object
    if(ncol(x$adata$data) == 0)
        x$adata$data = data.frame(x = rep(' ', nrow(x$adata$data)))
    if(!is.null(x$mlfit))
        p = ml.plot(x$mlfit, x$code, x$model, x$adata, x$sdata)
    if(!is.null(x$stanfit)){
        vs = names(x$adata$data)
        for(rpar in x$random)
            for(rg in rpar)
                vs = setdiff(vs, rg$group.name)
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
            p = ggplot(dfrocs, aes(cumfr, stim2$cumfr)) +
                geom_line(aes(x = cumfr.fit, y = stim2$cumfr.fit), lty = 2) +
                geom_errorbar(aes(ymin = stim2$cumfr.lo, ymax = stim2$cumfr.hi, x = cumfr.fit), width = 0.02) +
                geom_errorbarh(aes(xmin = cumfr.lo, xmax = cumfr.hi, y = stim2$cumfr.fit), height = 0.02) +
                geom_point() +
                labs(x = 'p(F)', y = 'p(H)') +
                coord_fixed() +
                facet_wrap(~f)
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
                p = ggplot(dfa, aes(response, count / n), group = stimulus) +
                    geom_errorbar(aes(ymin = count.lo / n, ymax = count.hi / n, lty = stimulus), width = 0.2) +
                    geom_line(aes(y = count.fit / n, lty = stimulus)) +
                    geom_point(aes(pch = stimulus)) +
                    labs(x = 'Response', y = 'Frequency', pch = 'Stimulus', lty = 'Stimulus') +
                    facet_wrap(~f) + theme_minimalist()
            }else{
                p = ggplot(dfa, aes(response, count / n, color = stimulus), group = stimulus) +
                    geom_errorbar(aes(ymin = count.lo / n, ymax = count.hi / n), width = 0.2) +
                    geom_line(aes(y = count.fit / n)) +
                    geom_point(aes(pch = stimulus)) +
                    labs(x = 'Response', y = 'Frequency', color = 'Stimulus', pch = 'Stimulus') +
                    facet_wrap(~f)
            }
            p
        }
    }
    p
}

ml.plot = function(mlfit, model_code, model, adata, sdata){
    mp = matrix(mlfit$par[grep('multinomial_p', names(mlfit$par))], nrow = max(1, nrow(adata$data)))
    obs = t(apply(sdata$counts, 1, function(x)(x / sum(x))))
    if(model %in% c('sdt', 'uvsdt', 'metad'))
        stim = as.factor(adata$stimulus)
    f = as.factor(adata$data[,1])
    if(ncol(adata$data) > 1)
        for(i in 2:ncol(adata$data))
            f = f:as.factor(adata$data[,i])
    title = sprintf('Lines = model: %s, points = observed response distribution', model)
    yl = 'p(Response)'
    if(model %in% c('sdt', 'uvsdt', 'metad')){
        xl = 'Response = Decision + Rating'
        p = ggplot(data.frame(f = rep(f, ncol(mp)), stim = rep(stim, ncol(mp)), th = rep(1:ncol(mp), each = length(f)),
                              mp = as.vector(mp), obs = as.vector(obs)),
                   aes(th, obs, group = stim, color = stim)) +
            geom_line(aes(y = mp)) +
            geom_point() +
            labs(title = title, color = 'Stimulus') + xlab(xl) + ylab(yl) +
            facet_wrap(~ f)
    }else{
        xl = 'Response = Rating'
        p = ggplot(data.frame(f = rep(f, ncol(mp)), th = rep(1:ncol(mp), each = length(f)),
                          mp = as.vector(mp), obs = as.vector(obs)),
                   aes(th, obs)) +
            geom_line(aes(y = mp)) +
            geom_point() +
            labs(title = title) +
            xlab(xl) +
            ylab(yl) +
            facet_wrap(~ f)
    }
    p
}

