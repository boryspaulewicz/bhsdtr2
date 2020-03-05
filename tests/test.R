## devtools::install('~/cs/code/r/bhsdtr2')
library(bhsdtr2)
library(rstan)
gabor$r = with(gabor, combined.response(stim, rating, acc))
gabor$r2 = with(gabor, combined.response(stim, accuracy = acc))

(m0 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r2 ~ stim,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),],
             method = F))
m0 = fit(m0)
cat(m0$code, file = '~/cs/code/r/bhsdtr2/tests/m0code.stan')

(m1.1 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r2 ~ stim,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),]))
## m1.1 = fit(m1.1, sample.prior = T)
samples(m1.1, 'dprim')
samples(m1.1, 'thr')

(m1.2 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),]))
(res = samples(m1.2, 'dprim'))
(res = samples(m1.2, 'thr'))
plot(m1.2)

## Sprawdzamy samples na poziomie indywidualnym
(m2 = bhsdtr(c(dprim ~ duration + (duration | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             method = 'stan'))
(res = samples(m2, 'dprim'))
(res = samples(m2, 'thr'))

dprim = samples(m2, 'dprim', group = 1)
thr = samples(m2, 'thr', group = 1, include.vars = 'duration')
all(attr(dprim, 'data') == attr(thr, 'data'))
df = attr(dprim, 'data')
df$dprim = apply(dprim[,1,], 2, mean)
df$thr = apply(thr[,4,], 2, mean)
## T³umaczenie miêdzy indeksem grupy a oryginaln± nazw± poziomu grupy
gr = unique(data.frame(index = m2$sdata$delta_group_1, original = m2$sdata$delta_group_1_original))
rownames(gr) = as.character(gr$index)
df$id.orig = gr[as.character(df$id), 'original']
df$acc.fit = acc.to.delta(delta = log(df$dprim), crit = df$thr)
df$acc = NA
for(i in 1:nrow(df))
    df$acc[i] = mean(gabor[(gabor$order == 'DECISION-RATING') & gabor$id == df$id.orig[i] & gabor$duration == df$duration[i], 'acc'])
plot(acc ~ acc.fit, df, main = sprintf('%.2f', cor(df$acc, df$acc.fit)))
## r = .94

## merging chains
(m3 = bhsdtr(c(dprim ~ order + (1 | id), thr ~ order + (1 | id)), r ~ stim,
             gabor[(gabor$duration == '32 ms'),],
             method = 'stan', chains = 2, iter = 2000, warmup = 1000))

merged.extract = function(m, par, group = NULL){x
    if(is.null(group)){
        par.name = sprintf('%s_fixed', par)
    }else{
        par.name = sprintf('%s_random_%d', par, group)
    }
    size = m$sdata[[sprintf('%s_size', par)]]
    res1 = extract(m$stanfit, par.name, permute = F)
    array(res1, c(dim(res1)[1] * dim(res1)[2], size, dim(res1)[3] / size))
}

res1 = merged.extract(m3, 'gamma', group = 1)
res2 = extract(m3$stanfit, 'gamma_random_1', permute = F)
iter = dim(res2)[1]
size = dim(res1)[2]
for(chain in 1:(dim(res2)[2])){
    for(cnd in 1:(dim(res1)[3])){
        for(d in 1:size)
            print(table(res1[1:iter + (chain - 1) * iter, d, cnd] == res2[,chain, d + (cnd - 1) * size]))
    }
}

plot(m2)

## Test set.prior
(m3 = bhsdtr(c(dprim ~ duration + (1 | id), thr ~ duration + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             method = F))

priors = list(delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))
m3.p = set.prior(m3, delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))
all(m3.p$sdata$delta_prior_fixed_mu == priors$delta$mu) &
all(m3.p$sdata$delta_prior_fixed_sd == priors$delta$sd) &
all(m3.p$sdata$delta_prior_scale_1 ==  priors$delta$scale[['1']]) &
all(m3.p$sdata$delta_prior_nu_1 ==  priors$delta$nu[['1']])

m3.p = set.prior(m3, delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))

(m3 = bhsdtr(c(dprim ~ duration + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             fit_method = 'stan'))
res1 = samples(m3, 'dprim')

(m4 = bhsdtr(c(mean ~ 1, thr ~ stim + (stim | id)), r ~ 1,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),],
             fit_method = 'stan'))
samples(m4, 'thr')
summary((x = samples(m4, 'thr')), digits = 3)
plot(m4)

######################################################################
## Testy

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

fit = 'ml'
models = c('sdt', 'uvsdt', 'metad')
links = c('softmax', 'log_distance', 'log_ratio', 'parsimonious', 'twoparameter')

test_fixed_models = function(method = 'jmap', links = c('softmax', 'log_distance', 'log_ratio', 'parsimonious', 'twoparameter'),
                            models = c('sdt', 'uvsdt')){ ## 'metad'
    d = sim_sdt(n = 10000)
    ## sd_ratio = 2 checked -> theta ~= log(2)
    table(d[, c('stim', 'r')])
    res = expand.grid(par = c('dprim', paste('c', 1:length(attr(d, 'criteria')), sep = '')),
                      link = links, model = models, true = NA, est = NA)
    for(model in levels(res$model)){
        for(link in levels(res$link)){
            print(sprintf('Fitting model %s with link %s', model, link))
            if(model == 'uvsdt'){
                m = bhsdtr(c(dprim ~ 1, thr ~ 1, sdratio ~ 1), r ~ stim, d, list(gamma = link), method = method)
            }else if(model == 'sdt'){
                m = bhsdtr(c(dprim ~ 1, thr ~ 1), r ~ stim, d, list(gamma = link), method = method)
            }else if(model == 'metad'){
                m = bhsdtr(c(metad ~ 1, thr ~ 1), r ~ stim, d, list(gamma = link), method = method)
            }
            res[(res$link == link) & (res$model == model), 'est'] = c(apply(samples(m, 'dprim')[,1,, drop = F], 3, mean),
                                                                      apply(samples(m, 'thr'), 2, mean))
            res[(res$link == link) & (res$model == model), 'true'] = c(dprim = attr(d, 'dprim'), attr(d, 'criteria'))
            rm(m)
            gc()
        }
    }
    title = sprintf('cor est and true: %f, method: %s', cor(res$est, res$true, use = 'pairwise.complete.obs'), method)
    print(title)
    res$par_type = c('criteria', 'dprim')[(as.character(res$par) == 'dprim') + 1]
    print(ggplot(res, aes(true, est, group = par_type, color = par_type)) +
          geom_abline(slope = 1, intercept = 0) +
          geom_point(size = 5) +
          facet_grid(model ~ link) +
          labs(title = title))
    res
}

test_fixed_models(models = 'sdt', links = 'parsimonious')
test_fixed_models('stan', models = 'sdt', links = 'parsimonious')
test_fixed_models(models = 'metad', links = 'parsimonious')
test_fixed_models(models = 'sdt')
test_fixed_models(models = 'uvsdt')
test_fixed_models(models = 'metad')

(m1 = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (order | id)),
             r ~ stim,
             gabor))
plot(m1)

(m2 = bhsdtr(c(dprim ~ -1 + duration:order + (-1 + duration | id), thr ~ order + (order | id)),
             r ~ stim,
             gabor, fit = 'stan'))
plot(m2)

## bhsdtr 1
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'order', 'id'))
sdata = make_stan_data(adata, list(delta = ~ -1 + duration:order, gamma = ~ order), list(list(delta = ~ -1 + duration, gamma = ~ order, group = ~ id)),
                       gamma_link = 'log_distance')
model_code = make_stan_model(list(list(delta = ~ -1 + duration, gamma = ~ order, group = ~ id)), links$gamma)
res1 = optimizing(stan_model(model_code = model_code), sdata, init = 0)
if(res1$return_code)print('not ok!')
