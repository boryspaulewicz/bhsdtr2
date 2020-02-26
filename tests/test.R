library(bhsdtr2)
library(rstan)

gabor$r = with(gabor, combined.response(stim, rating, acc))

(m1 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),]))

(res = samples(m1, 'dprim'))
plot(m1)

(m2 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING') & (gabor$duration == '32 ms'),],
             fit_method = 'stan'))
(res = samples(m2, 'dprim'))
(res = samples(m2, 'thr'))
plot(m2)

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

######################################################################


samples(m, 'dprim')

## Testy
fit = 'ml'
models = c('sdt', 'uvsdt', 'metad')
links = c('softmax', 'log_distance', 'log_ratio', 'parsimonious', 'twoparameter')

test_fixed_models = function(fit = 'ml', links = c('softmax', 'log_distance', 'log_ratio', 'parsimonious', 'twoparameter'),
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
                m = bhsdtr(c(dprim ~ 1, thr ~ 1, sdratio ~ 1), r ~ stim, d, list(gamma = link), fit = fit)
            }else if(model == 'sdt'){
                m = bhsdtr(c(dprim ~ 1, thr ~ 1), r ~ stim, d, list(gamma = link), fit = fit)
            }else if(model == 'metad'){
                m = bhsdtr(c(metad ~ 1, thr ~ 1), r ~ stim, d, list(gamma = link), fit = fit)
            }
            res[(res$link == link) & (res$model == model), 'est'] = c(apply(t(samples(m, 'dprim')[,1,]), 2, mean), apply(samples(m, 'thr'), 2, mean))
            res[(res$link == link) & (res$model == model), 'true'] = c(dprim = attr(d, 'dprim'), attr(d, 'criteria'))
            rm(m)
            gc()
        }
    }
    print(sprintf('cor est and true: %f', cor(res$est, res$true, use = 'pairwise.complete.obs')))
    res$par_type = c('criteria', 'dprim')[(as.character(res$par) == 'dprim') + 1]
    print(ggplot(res, aes(est, true, group = par_type, color = par_type)) +
          geom_abline(slope = 1, intercept = 0) +
          geom_point(size = 5) +
          facet_grid(model ~ link))
    res
}

test_fixed_models(models = 'sdt')
test_fixed_models(models = 'uvsdt')
test_fixed_models(models = 'metad', links = 'parsimonious')
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
