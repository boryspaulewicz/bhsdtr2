## devtools::install('~/Dropbox/CS/code/r/bhsdtr2')
## devtools::document('~/cs/code/r/bhsdtr2')
source('~/Dropbox/CS/code/r/bhsdtr2/tests/utils.R')
library(bhsdtr2)
## library(bhsdtr)
library(rstan)
library(bridgesampling)
gabor$r = with(gabor, combined.response(stim, rating, acc))
gabor$r2 = with(gabor, combined.response(stim, accuracy = acc))

######################################################################
## plot jmap

(m = bhsdtr(c(dprim ~ duration + (duration | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[(gabor$order == 'DECISION-RATING'),],
            prior = list(gamma_prior_random_scale_1 = 0)))

(m = bhsdtr(c(dprim ~ -1 + duration + (-1 + duration | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[(gabor$order == 'DECISION-RATING'),],
            links = list(delta = 'id_log')))
## Najs�absze korelacje z pozosta�ymi efektami losowymi ma odchylenie od �redniego kryterium
## �rodkowego.
round(apply(matrix(res[grepl('^corr_gamma', names(res))], nrow = 7), 2, mean), 2)
## 0.21 0.26 0.24 0.13 0.22 0.24 0.23

## Korelacje wykazuj� efekt symetrii - odchylenia dla 1go progu koreluj� najsilniej z odchyleniami
## dla 7go progu, te dla 2go najsilniej z tymi dla 6go, i tak dalej: i najsilniej z 8 - i. To ma
## sens, bo skrajne oceny to w pewnym sensie te same oceny (ten sam rating, ale r�ne decyzje)
apply(matrix(res[grepl('^corr_gamma', names(res))], nrow = 7), 2, function(x)which(x == max(x[x != 1])))
## 7 6 5 3 5 2 1
plot(m)
plot(m, vs = 'duration')

samples(m, 'dprim')
samples(m, 'thr')

(m = bhsdtr(c(dprim ~ -1 + duration + (-1 + duration | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[(gabor$order == 'DECISION-RATING'),],
            links = list(gamma = 'id_log'),
            method = ''))
cat(m$code, file = 'gamma_id_log.stan')

######################################################################
## simple hierarchical sdt model

## bhsdtr1 fit
adata = aggregate_responses(gabor[(gabor$order == 'DECISION-RATING' & gabor$duration == '32 ms'),], 'stim', 'r', c('id'))
sdata = make_stan_data(adata, list(delta = ~ 1, gamma = ~ 1), list(list(delta = ~ 1, gamma = ~ 1, group = ~ id)),
                       gamma_link = 'log_distance')
model_code = make_stan_model(list(list(delta = ~ 1, gamma = ~ 1, group = ~ id)), gamma_link = 'log_distance')
res1 = optimizing(stan_model(model_code = model_code), sdata, init = 0)
if(res1$return_code)print('not ok!')
round(exp(res1$par[['delta_fixed[1,1]']]), 2)
## .91

(m = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[(gabor$order == 'DECISION-RATING' & gabor$duration == '32 ms'),]))
samples(m, 'dprim')
## samples: 1, estimates rounded to 2 decimal places
##  dprim.1
##     1.18
samples(m, 'thr')
## samples: 1, estimates rounded to 2 decimal places
##  thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
##  -2.00 -1.33 -0.77  0.08  0.62  1.07  1.87
plot(m)

## bhsdtr1 with bhsdtr2 default priors
sdata2 = sdata
sdata2$delta_fixed_mu = m$sdata$delta_prior_fixed_mu
sdata2$delta_fixed_sd = m$sdata$delta_prior_fixed_sd
sdata2$delta_sd_scale_1 = m$sdata$delta_prior_random_scale_1
sdata2$gamma_fixed_mu = m$sdata$gamma_prior_fixed_mu
sdata2$gamma_fixed_sd = m$sdata$gamma_prior_fixed_sd
sdata2$gamma_sd_scale_1 = m$sdata$gamma_prior_random_scale_1
res2 = optimizing(stan_model(model_code = model_code), sdata2, init = 0)
if(res2$return_code)print('not ok!')
round(exp(res2$par[['delta_fixed[1,1]']]), 2)
## 1.18, which is the same

(m = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[(gabor$order == 'DECISION-RATING' & gabor$duration == '32 ms'),],
            method = 'stan'))
samples(m, 'dprim')
## samples: 21000, estimates rounded to 2 decimal places
##  dprim.1
##     0.93
samples(m, 'thr')
## samples: 21000, estimates rounded to 2 decimal places
##  thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
##  -1.88 -1.46 -0.93  0.08  0.67  1.08  1.78

######################################################################
## bf bridgesampling

(m0 = bhsdtr(c(dprim ~ -1 + duration + (-1 + duration | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             links = list(delta = 'id_log'),
             method = 'stan'))

(m1 = bhsdtr(c(dprim ~ -1 + duration + (-1 + duration | id), thr ~ -1 + duration + (-1 + duration | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             links = list(delta = 'id_log'),
             method = 'stan'))

summary((b0 = bridge_sampler(m0$stanfit)))
## std links
##
## perc. err 49%
##
## id_log
##
## perc. err  77%

summary((b0 = bridge_sampler(m0$stanfit, method = 'warp3', repetitions = 10, cores = 7, silent = T)))
summary((b1 = bridge_sampler(m1$stanfit, method = 'warp3', repetitions = 10, cores = 7, silent = T)))
bf(b0, b1)
## std links
##
## Estimated Bayes factor (based on medians of log marginal likelihood estimates)
##  in favor of b0 over b1: 358516.24015
## Range of estimates: 4219.02809 to 3400546.54379
## Interquartile range: 1622755.75089
##
## id_log link
##
## Estimated Bayes factor (based on medians of log marginal likelihood estimates)
##  in favor of b0 over b1: 60813763255323408.00000
## Range of estimates: 94527492123870.76562 to 17544479217032108032.00000
## Interquartile range: 909534897459633664.00000

######################################################################
## cumulative

m = cumulative(r ~ stim + (stim | id),
               gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',])
plot(m)
samples(m, 'thr')

######################################################################
## Testing participant-level estimates from the samples function

(m2.jmap = bhsdtr(c(dprim ~ duration + (duration | id), thr ~ 1 + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),]))
samples(m2.jmap, 'dprim')
sdt.acc.plot(m2.jmap, 1)

######################################################################
## set.prior

(m3 = bhsdtr(c(dprim ~ duration + (1 | id), thr ~ duration + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             method = F))

(m3.p = bhsdtr(c(dprim ~ duration + (1 | id), thr ~ duration + (1 | id)), r ~ stim,
             gabor[(gabor$order == 'DECISION-RATING'),],
             method = F, prior = list(delta = list(mu = 123, sd = 321))))
m3.p$sdata$delta_prior_fixed_mu
m3.p$sdata$delta_prior_fixed_sd

priors = list(delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))
m3.p = set.prior(m3, delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))
all(m3.p$sdata$delta_prior_fixed_mu == priors$delta$mu) &
all(m3.p$sdata$delta_prior_fixed_sd == priors$delta$sd) &
all(m3.p$sdata$delta_prior_scale_1 ==  priors$delta$scale[['1']]) &
all(m3.p$sdata$delta_prior_nu_1 ==  priors$delta$nu[['1']])

m3.p = set.prior(m3, delta = list(mu = 123, sd = 234, scale = list('1' = 321), nu = list('1' = 333)))

######################################################################
## models and links tests

## Test meta-d' estimation
n = 100000
dprim = 2
metad = 1
criteria = c(-2.1, -1.4, -.7, 0, .7, 1.4, 2.1)

cum1.s1 = pnorm(c(criteria, +Inf) + dprim / 2)
probs1.s1 = cum1.s1 - c(0, cum1.s1[-length(cum1.s1)])
cum1.s2 = pnorm(c(criteria, +Inf) - dprim / 2)
probs1.s2 = cum1.s2 - c(0, cum1.s2[-length(cum1.s2)])

cum2.s1 = pnorm(c(criteria, +Inf) + metad / 2)
probs2.s1 = cum2.s1 - c(0, cum2.s1[-length(cum2.s1)])
cum2.s2 = pnorm(c(criteria, +Inf) - metad / 2)
probs2.s2 = cum2.s2 - c(0, cum2.s2[-length(cum2.s2)])
K = length(criteria) + 1
Kb2 =  K / 2
hit1 = sum(probs1.s2[(Kb2 + 1):K])
fa1 = sum(probs1.s1[(Kb2 + 1):K])
hit2 = sum(probs2.s2[(Kb2 + 1):K])
fa2 = sum(probs2.s1[(Kb2 + 1):K])
## Normalization
probs2.s2[(Kb2 + 1):K] = probs2.s2[(Kb2 + 1):K] * hit1 / hit2
probs2.s1[(Kb2 + 1):K] = probs2.s1[(Kb2 + 1):K] * fa1 / fa2
probs2.s2[1:Kb2] = probs2.s2[1:Kb2] * (1 - hit1) / (1 - hit2)
probs2.s1[1:Kb2] = probs2.s1[1:Kb2] * (1 - fa1) / (1 - fa2)

d = data.frame(stim = rep(1:2, each = n))
d$r[d$stim == 1] = apply(rmultinom(n, 1, probs1.s1), 2, function(x)which(x == 1))
d$r[d$stim == 2] = apply(rmultinom(n, 1, probs1.s2), 2, function(x)which(x == 1))
d$r2[d$stim == 1] = apply(rmultinom(n, 1, probs2.s1), 2, function(x)which(x == 1))
d$r2[d$stim == 2] = apply(rmultinom(n, 1, probs2.s2), 2, function(x)which(x == 1))

m1 = bhsdtr(c(dprim ~ 1, thr ~ 1), r ~ stim, d)
samples(m1, 'dprim')
## Ok
m2 = bhsdtr(c(metad ~ 1, thr ~ 1), r2 ~ stim, d)
samples(m2, 'dprim')
## Ok

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
## test_fixed_models('stan', models = 'sdt', links = 'parsimonious')
test_fixed_models(models = 'metad', links = 'parsimonious')
test_fixed_models(models = 'sdt')
test_fixed_models(models = 'uvsdt')
test_fixed_models(models = 'metad')

######################################################################
## bhsdtr 1

adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'order', 'id'))
sdata = make_stan_data(adata, list(delta = ~ -1 + duration:order, gamma = ~ order), list(list(delta = ~ -1 + duration, gamma = ~ order, group = ~ id)),
                       gamma_link = 'log_distance')
model_code = make_stan_model(list(list(delta = ~ -1 + duration, gamma = ~ order, group = ~ id)), links$gamma)
res1 = optimizing(stan_model(model_code = model_code), sdata, init = 0)
if(res1$return_code)print('not ok!')
