######################################################################
## This is the code that was used to perform all the analyses and
## produce all the plots in https://psyarxiv.com/e7a3x/

## WARNING You will probably have to alter the model.path function
## which is defined in this script - this function is appropriate to
## my setup but it will almost surely not be appropriate for yours.

library(rstanarm)
options(mc.cores = parallel::detectCores())
library(bhsdtr2)
library(lme4)
library(ggplot2)
library(ggthemes)
library(coda)
data(gabor)
gabor$r = combined.response(gabor$stim, gabor$rating, gabor$acc)
## We will be using this subset of the gabor dataset
d = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING', ]
## We have to make sure that the results are reproducible
seed = 12345

my_theme = function()theme_few() ## geom_rangeframe() + theme_tufte(ticks = F)

gamma_link = 'log_distance'

## The fitted models are saved for later
model.path = function(name)sprintf('E:/temp/cumulative_%s/%s.rds', gamma_link, name)
read.model = function(name)readRDS(model.path(name))
model.is.saved = function(name)file.exists(model.path(name))
save.model = function(model){
    name = deparse(substitute(model))
    saveRDS(model, file = model.path(name))
}
fun = function(x)quantile(x, c(.025, .975))
funhpd90 = function(x)HPDinterval(as.mcmc(x), .90)
funhpd95 = function(x)HPDinterval(as.mcmc(x), .95)
is.in = function(s, true){
    fun = function(x)quantile(x, c(.025, .975))
    qs = t(apply(s, 3, fun))
    qs[,1] <= true & qs[,2] >= true
}

## rt stands for random thresholds
if(!model.is.saved('m.sdt.rt')){
    m.sdt.rt = bhsdtr(c(dprim ~ 1 + (1|id), thr ~ 1 + (1|id)), r ~ stim, d, links = list(gamma = gamma_link), seed = seed)
    ## plot(m.sdt.rt, vs = 'id', bw = T) + my_theme()
    ## ggsave('m_sdt_rt.pdf', width = 16, height = 7)
    ## plot(m.sdt.rt, vs = 'id', bw = T, type = 'roc') + my_theme()
    ## ggsave('m_sdt_rt_roc.pdf', width = 7, height = 7)
    m.sdt.rt = fit(m.sdt.rt, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.sdt.rt)
}else{
    m.sdt.rt = read.model('m.sdt.rt')
}

## The general cumulative model
if(!model.is.saved('m.cm')){
    m.cm = cumulative(r ~ stim + (stim | id), d, method = 'stan', links = list(gamma = gamma_link), seed = seed,
                      chains = 7, iter = 5000, warmup = 2000)
    save.model(m.cm)
}else{
    m.cm = read.model('m.cm')
}

######################################################################
## All the following models are fitted to simulated data

## Counts simulated from the known realistic model
dprim = samples(m.sdt.rt, 'dprim', 1, method = 'jmap')[1,,]
thr = samples(m.sdt.rt, 'thr', 1, method = 'jmap')[1,,]
sim.true.sdt = function(dprim, thr, n){
    set.seed(seed)
    ds = NULL
    for(i in 1:(nrow(m.sdt.rt$sdata$counts) / 2)){
        res = sim.sdt(n, dprim[i], thr[,i])
        res$id = i
        if(is.null(ds)){
            ds = res
        }else{
            ds = rbind(ds, res)
        }
    }
    ds
}
ds = sim.true.sdt(dprim, thr, 100)
ds48 = sim.true.sdt(dprim, thr, 24)

ds$resp = as.numeric(ds$r > 4) + 1
ds$dprim = dprim[ds$id]
ds$acc = as.numeric(ds$stim == ds$resp)
res = aggregate(acc ~ id + dprim, ds, sum)
cor(res[,c('dprim', 'acc')])
## .9

## The true model
if(!model.is.saved('m.true')){
    m.true = bhsdtr(c(dprim ~ 1 + (1|id), thr ~ 1 + (1|id)), r ~ stim, ds, links = list(gamma = gamma_link), seed = seed)
    m.true = fit(m.true, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.true)
    ## plot(m.true, fit = 'jmap')
    ## ggsave('m.true.jpg')
}else{
    m.true = read.model('m.true')
}

## The true cumulative model fitted to data for stimulus class 1.
if(!model.is.saved('m.half.1')){
    m.half.1 = cumulative(r ~ 1 + (1 | id), ds[ds$stim == 1,], links = list(gamma = gamma_link), seed = seed)
    m.half.1 = fit(m.half.1, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.half.1)
}else{
    m.half.1 = read.model('m.half.1')
}

if(!model.is.saved('m.half.2')){
    m.half.2 = cumulative(r ~ 1 + (1 | id), ds[ds$stim == 2,], links = list(gamma = gamma_link), seed = seed)
    m.half.2 = fit(m.half.2, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.half.2)
}else{
    m.half.2 = read.model('m.half.2')
}

## The fixed thresholds (ft) model
if(!model.is.saved('m.sdt.ft')){
    m.sdt.ft = bhsdtr(c(dprim ~ 1 + (1|id), thr ~ 1), r ~ stim, ds, links = list(gamma = gamma_link), seed = seed)
    m.sdt.ft = fit(m.sdt.ft, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.sdt.ft)
    ## plot(m.sdt.ft, bw = T) + my_theme()
    ## ggsave('m_sdt_ft.pdf', width = 7, height = 7)
    ## plot(m.sdt.ft, bw = T, type = 'roc') + my_theme()
    ## ggsave('m_sdt_ft_roc.pdf', width = 7, height = 7)
}else{
    m.sdt.ft = read.model('m.sdt.ft')
}

## The fixed threshold (ft) model without ratings
if(!model.is.saved('m.sdt.c')){
    ds$r2 = as.numeric(ds$r > 4)
    m.sdt.c = bhsdtr(c(dprim ~ 1 + (1|id), thr ~ 1), r2 ~ stim, ds, links = list(gamma = gamma_link), seed = seed)
    m.sdt.c = fit(m.sdt.c, method = 'stan', seed = seed, chains = 7, iter = 5000, warmup = 2000)
    save.model(m.sdt.c)
}else{
    m.sdt.c = read.model('m.sdt.c')
}

######################################################################
## Analyses

s.sdt.rt = samples(m.sdt.rt, 'dprim', 1)
s.true = samples(m.true, 'dprim', 1)
s.sdt.ft = samples(m.sdt.ft, 'dprim', 1)
s.sdt.c = samples(m.sdt.c, 'dprim', 1)
s.half.1 = samples(m.half.1, 'thr', 1)

cor(apply(s.true, 3, mean), dprim)
## 0.9798475
cor(apply(s.sdt.rt, 3, mean), dprim)
## 0.9899008
cor.test(apply(s.sdt.c, 3, mean), dprim)
## .96
cor.test(apply(s.sdt.ft, 3, mean), dprim)
## 0.777477, t(27) = 6.42, p < .001, CI = .57-.89, softmax: 0.7752282
m = lmer(r ~ stim + (stim | id), d)
cor(ranef(m)$id[,2], dprim)
## 0.7837598

var(apply(s.true, 3, mean))
## 0.942467
var(apply(s.sdt.ft, 3, mean))
## 0.6985466
var(apply(s.true, 3, mean)) / var(apply(s.sdt.ft, 3, mean))
## 1.349183
var(apply(s.sdt.c, 3, mean))
## 0.7737335
var(apply(s.true, 3, mean)) / var(apply(s.sdt.c, 3, mean))
## 1.218077

sd(samples(m.sdt.ft, 'dprim')[,1,1]) / sd(samples(m.true, 'dprim')[,1,1])
1 - mean(apply(samples(m.sdt.ft, 'thr')[,,1], 2, sd) / apply(samples(m.true, 'thr')[,,1], 2, sd))
round(mean(samples(m.sdt.rt, 'dprim')[,1,1]), 2)
## .94
round(fun(samples(m.true, 'dprim')[,,1]), 2)
## 2.5% 97.5% 
## 0.62  1.38 
round(fun(samples(m.sdt.ft, 'dprim')[,,1]), 2)
## 2.5% 97.5%
## 0.50 1.07

round(apply(samples(m.sdt.rt, 'thr')[,,1], 2, mean), 2)
##       thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7 
##       -1.88 -1.46 -0.92  0.08  0.67  1.08  1.78 
round(apply(samples(m.true, 'thr')[,,1], 2, fun), 2)
##       thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
## 2.5%  -2.48 -1.88 -1.26 -0.07  0.46  0.84  1.57
## 97.5% -1.55 -1.17 -0.67  0.28  0.89  1.35  2.26
round(apply(samples(m.sdt.ft, 'thr')[,,1], 2, fun), 2)
##       thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
## 2.5%  -1.87 -1.42 -0.97  0.05  0.69  1.09  1.75
## 97.5% -1.75 -1.33 -0.89  0.12  0.77  1.18  1.86
##                *                *     *

## A simple measure of estimate bias: distance from the point estimate
## to the true value in units of standard deviation of posterior
## samples
dp = data.frame(y = c((apply(s.true, 3, mean) - dprim) / apply(s.true, 3, sd),
(apply(s.sdt.ft, 3, mean) - dprim) / apply(s.sdt.ft, 3, sd)))
dp$model = rep(c('true', 'simplified'), each = dim(s.true)[3])
dp$i = 1:dim(s.true)[3]
ggplot(dp, aes(i, y, group = model, shape = model, color = model)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 0)
## ggsave('bias.dprim.jpg')

## A plot of point and interval estimates vs true values
dp = data.frame(dprim = c(dprim, dprim),
                pointest = c(apply(s.true, 3, mean),
                             apply(s.sdt.ft, 3, mean)),
                lo = c(apply(s.true, 3, function(x)quantile(x, .025)),
                       apply(s.sdt.ft, 3, function(x)quantile(x, .025))),
                hi = c(apply(s.true, 3, function(x)quantile(x, .975)),
                       apply(s.sdt.ft, 3, function(x)quantile(x, .975))))
dp$model = rep(c('True model', 'Simplified model'), each = dim(s.true)[3])
dp$i = 1:dim(s.true)[3]

ggplot(dp, aes(i, pointest)) +
    ## geom_point() +
    geom_point(aes(y = dprim), size = 2) +
    geom_segment(aes(x = i, xend = i, y = pointest, yend = dprim), lty = 'dotted') + 
    geom_errorbar(aes(ymin = lo, ymax = hi)) +
    labs(y = 'Sensitivity (d\')', x = 'Participant index') +
    facet_grid(. ~ model) + theme_hc() + theme(strip.background = element_blank()) ## + my_theme() ## + theme(panel.grid.major = element_blank())

ggsave('dprim_point_interval_estimates.pdf', width = 16, height = 7)

mean(apply(s.true, 3, sd))
## 0.1744569
mean(apply(s.sdt.ft, 3, sd))
## 0.1361601
mean(apply(s.sdt.ft, 3, sd)) / mean(apply(s.true, 3, sd))
## 0.7804801 softmax: 0.7481904

dp = data.frame(y = c((apply(s.true, 3, mean) - dprim) / apply(s.true, 3, sd),
(apply(s.sdt.ft, 3, mean) - dprim) / apply(s.sdt.ft, 3, sd)))
dp$model = rep(c('true', 'simplified'), each = dim(s.true)[3])
dp$i = 1:dim(s.true)[3]
## Average bias severity, 1 is fine
aggregate(y ~ model, dp, function(x)mean(abs(x)))
##        model        y
## 1 simplified 3.242385
## 2       true 0.895169

## SDT 1 thr
dp.c = data.frame(y = c((apply(s.true, 3, mean) - dprim) / apply(s.true, 3, sd),
(apply(s.sdt.c, 3, mean) - dprim) / apply(s.sdt.c, 3, sd)))
dp.c$model = rep(c('true', 'sdt.c'), each = dim(s.true)[3])
dp.c$i = 1:dim(s.true)[3]
## Average bias severity, 1 is fine
aggregate(y ~ model, dp.c, function(x)mean(abs(x)))
##   model         y
## 1 sdt.c 0.8984822
## 2  true 0.8951690

## softmax
##        model         y
## 1 simplified 3.2744274
## 2       true 0.8784176

## CI coverage
aggregate(y ~ model, dp, function(x)mean(abs(x) > 1.96))
## 1 simplified 0.55172414
## 2       true 0.03448276
aggregate(y ~ model, dp.c, function(x)mean(abs(x) > 1.96))
##   model          y
## 1 sdt.c 0.10344828
## 2  true 0.03448276

## softmax
##        model          y
## 1 simplified 0.55172414
## 2       true 0.03448276

## CI coverage
mean(is.in(s.sdt.rt, dprim))
## 0.8965517 softmax: 0.9310345
mean(is.in(s.true, dprim))
## 0.9310345 softmax: 0.9310345
1 - mean(is.in(s.true, dprim))
mean(is.in(s.sdt.ft, dprim))
## 0.4482759 softmax: 0.4482759
1 - mean(is.in(s.sdt.ft, dprim))
mean(is.in(s.sdt.c, dprim))
## 0.8965517

## ## Group-level fit
## plot(m.true, vs = 1)
## ggsave('m.true.fit.jpg')
## plot(m.sdt.ft, vs = 1)
## ggsave('m.sdt.ft.fit.jpg')
## plot(m.true)
## ggsave('m.true.id.fit.jpg')
## plot(m.sdt.ft)
## ggsave('m.sdt.ft.id.fit.jpg')

######################################################################
## Testing causal assumptions

s.cm = samples(m.cm, 'thr', 1)
s.cm.48 = samples(m.cm.48, 'thr', 1)

## The effect of the stimulus class on the middle threshold is very
## close to dprim, although the effect on the average threshold
## position is not that close. That may be because the estimates of
## the more outermost thresholds are less precise.
sel.test = function(s.cm, dprim){
    s0 = s.cm[,,grep(':0$', dimnames(s.cm)[[3]])]
    s1 = s.cm[,,grep(':1$', dimnames(s.cm)[[3]])]
    dist.4 = dist.mean = matrix(nrow = dim(s0)[1], ncol = dim(s0)[3])
    for(i in 1:(dim(s0)[3])){
        dist.4[,i] = s0[,4,i] - s1[,4,i]
        dist.mean[,i] = apply(s0[,,i] - s1[,,i], 1, mean)
    }
    cor(apply(dist.4, 2, mean), apply(dist.mean, 2, mean))
    ## 0.3504015
    ## lmrob(apply(dist.4, 2, mean) ~ apply(dist.mean, 2, mean))
    ## ## Coefficients:
    ## ##               (Intercept)  apply(dist.mean, 2, mean)  
    ## ##                    0.2154                     0.8729  
    lm(apply(dist.4, 2, mean) ~ apply(dist.mean, 2, mean))
    ## Coefficients:
    ##               (Intercept)  apply(dist.mean, 2, mean)  
    ##                    1.0843                     0.1118  
    plot(apply(dist.4, 2, mean) ~ apply(dist.mean, 2, mean), xlim = c(0,3), ylim = c(0,3))
    cor(dprim, apply(dist.4, 2, mean))
    ## 0.9438046
    cor(dprim, apply(dist.mean, 2, mean))
    ## 0.370617
    mean(t(apply(dist.4, 2, funhpd95))[,1] <= dprim & t(apply(dist.4, 2, funhpd95))[,2] >= dprim)
    ## 0.9310345 softmax: 0.9310345
    mean(t(apply(dist.mean, 2, funhpd95))[,1] <= dprim & t(apply(dist.mean, 2, funhpd95))[,2] >= dprim)
    ## 0.9310345
    mean(apply(dist.4, 2, sd))
    ## 0.3390451 softmax: 0.7931034
    sc = s0
    for(i in 1:(dim(s0)[3])){
        ## sc[,,i] = s0[,,i] - s0[,4,i] - (s1[,,i] - s1[,4,i])
        sc[,,i] = s0[,,i] - apply(s0[,,i], 1, mean) - (s1[,,i] - apply(s1[,,i], 1, mean))
    }
    res = expand.grid(id = 1:(dim(s0)[3]), thr = 1:7, lo = NA, hi = NA)
    for(i in 1:nrow(res)){
        res[i, c('lo', 'hi')] = funhpd90(sc[, res$thr[i], res$id[i]])
    }
    with(res, mean(lo * hi > 0))
    ## 0.01970443
    ggplot(res, aes(as.factor(thr), group = id)) +
        geom_abline(intercept = 0, slope = 0, lty = 'dotted') + 
        geom_errorbar(aes(ymax = hi, ymin = lo, width = .2)) +
        labs(x = 'Threshold index', y = 'Difference between the distances from the average threshold position') +
        facet_wrap(~ id) + geom_rangeframe() + theme_tufte(ticks = F) ## my_theme()
}

sel.test(s.cm, dprim)
dimnames(s.cm.48)[[3]] = gsub(':1', ':0', dimnames(s.cm.48)[[3]])
dimnames(s.cm.48)[[3]] = gsub(':2', ':1', dimnames(s.cm.48)[[3]])
sel.test(s.cm.48, dprim)
ggsave('thresholds_contrast.pdf', width = 16, height = 7)
