## Próba rozdzielenia priorów i funkcji ³±cz±cych dla efektów
## ustalonych i losowych. Dodajemy nowy rodzaj funkcji ³±cz±cej dla
## parametrów dodatnich: id_log, który wymaga, ¿eby parametryzacja
## efektów ustalonych i losowych by³a separate.intercepts.

## Trzeba bêdzie porzuciæ dotychczasowego parsera i zacz±æ u¿ywaæ
## stringów wewn±trz funkcji tworz±cej kod modelu, albo dodaæ
## alternatywy dla obecno¶ci i nieobecno¶ci efektów losowych.

library(bhsdtr2)
library(rstan)
library(bayesplot)
library(ggplot2)

gabor$r = with(gabor, combined.response(stim, rating, acc))

m1 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)), r ~ stim,
            gabor[gabor$order == 'DECISION-RATING' & gabor$duration == '32 ms',])
samples(m1, 'dprim')
## dprim.1
##    0.91
samples(m1, 'thr')
samples: 1, estimates rounded to 2 decimal places
## thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
## -2.00 -1.33 -0.77  0.08  0.62  1.07  1.87
cat(m1$code, file = '~/cs/code/r/bhsdtr2/tests/model.stan')

## Dostosowujemy priory dla delta_fixed
m2 = m1
m2$code = paste(readLines('~/cs/code/r/bhsdtr2/tests/model2.stan'), collapse = '\n')
m2$sdata$delta_prior_fixed_mu[,1] = 1
m2$sdata$delta_prior_fixed_sd[,1] = 1
m2 = fit(m2, method = 'stan')
condition.specific.samples(m2, 'delta')
## delta.1
##    0.95 ## .92 gdy jmap
samples(m2, 'thr')
## samples: 1, estimates rounded to 2 decimal places
##  thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
##  -2.00 -1.33 -0.76  0.08  0.62  1.07  1.87
##
## stan
##
## samples: 21000, estimates rounded to 2 decimal places
##  thr.1 thr.2 thr.3 thr.4 thr.5 thr.6 thr.7
##  -1.89 -1.46 -0.93  0.08  0.67  1.09  1.78

round(rbind(condition.specific.samples(m1, 'gamma')[1,,1],
            apply(condition.specific.samples(m2, 'gamma')[,,1], 2, mean)), 2)
##      gamma.1 gamma.2 gamma.3 gamma.4 gamma.5 gamma.6 gamma.7
## [1,]   -0.41   -0.57   -0.16    0.08   -0.62    -0.8   -0.22
## [2,]   -0.41   -0.57   -0.17    0.08   -0.62    -0.8   -0.22
##
## stan
##
##      gamma.1 gamma.2 gamma.3 gamma.4 gamma.5 gamma.6 gamma.7
## [1,]   -0.41   -0.57   -0.16    0.08   -0.62   -0.80   -0.22
## [2,]   -0.89   -0.64    0.00    0.08   -0.54   -0.89   -0.37


## Wydaje siê dzia³aæ bardzo dobrze!

s = as.array(m2$stanfit)
mcmc_dens(s, pars = names(s[1,1,])[1:8])
## Piêknie
