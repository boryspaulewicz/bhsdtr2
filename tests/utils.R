## Nie wiem jak sobie poradzi z samym interceptem
sdt.acc.plot = function(model, group = NULL){
    if(!is.null(group)){
        group.name = model$random$delta[[group]]$group.name
    }else{
        group.name = NULL
    }
    vars = predictor.names(model, 'delta', 'gamma')
    dprim = samples(model, 'dprim', group = group, include = vars)
    thr = samples(model, 'thr', group = group, include = vars)
    df = samples.data(dprim)
    if(dim(dprim)[1] == 1){
        dprim = dprim[1,1,]
        crit = thr[1,model$sdata$Kb2,]
    }else{
        dprim = apply(dprim[,1,], 2, mean)
        crit = apply(thr[,model$sdata$Kb2,], 2, mean)
    }
    acc = paste(c(names(model$data), 'acc'), collapse = '_')
    model$data[[acc]] = as.numeric(((model$resp > model$sdata$Kb2) + 1) == model$stim)
    aggr.formula = as.formula(sprintf('%s ~ %s', acc, paste(c(vars, group.name), collapse = ' + ')))
    res = aggregate(aggr.formula, model$data, mean)
    res.rows = df.rows = ''
    for(v in c(vars, group.name)){
        res.rows = paste(res.rows, as.character(res[[v]]), sep = ':')
        df.rows = paste(df.rows, as.character(df[[v]]), sep = ':')
    }
    rownames(res) = res.rows
    acc.obs = paste(c(acc, 'obs'), collapse = '_')
    acc.fit = paste(c(acc, 'fit'), collapse = '_')
    df[, acc.obs] = res[df.rows, acc]
    df[, acc.fit] = acc.to.delta(delta = log(dprim), crit = crit)
    p = ggplot(df, aes_string(acc.fit, acc.obs)) +
        geom_point() + geom_abline(slope = 1, intercept = 0, lty = 2)
    if(!is.null(group))
        p = p + facet_wrap(as.formula(sprintf('~ %s', paste(vars, collapse = ' + '))))
    p = p + labs(title = sprintf('r = %.2f', cor(df[[acc.obs]], df[[acc.fit]])), x = 'Predicted accuracy', y = 'Observed accuracy')
    p
}

sim_sdt = function(n = 1, dprim = 1.5, criteria = c(-2.1, -1.4, -.7, 0, .7, 1.4, 2.1), sd_ratio = 1){
    which_bin = function(x, thr)min(which(x <= c(-Inf, thr, Inf)) - 1)
    d = data.frame(stim = rep(1:2, each = n), e = rnorm(n * 2), r = NA)
    d$e[d$stim == 2] = d$e[d$stim == 2] * sd_ratio
    d$e = d$e + .5 * dprim * c(-1, 1)[d$stim]
    for(i in 1:nrow(d))
        d$r[i] = which_bin(d$e[i], criteria)
    attr(d, 'dprim') = dprim
    attr(d, 'criteria') = criteria
    attr(d, 'sd_ratio') = sd_ratio
    d
}

