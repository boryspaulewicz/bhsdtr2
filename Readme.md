# This package is in the stage of active development

# NEWS 2020 02 17

[Here](ordinal/Readme.md) is a tutorial on hierarchical ordinal
polytomous models in bhsdtr, with an emphasis on SDT models (and a
rant about Item Response Theory).

# NEWS 2020 02 10

The bhsdtr package now handles hierarchical or non-hierarchical
versions of Selker, van den Bergh, Criss, and Wagenmakers's
parsimonious SDT model (see this
[paper](https://link.springer.com/article/10.3758/s13428-019-01231-3)).
It can be fitted by choosing the 'parsimonious' link function, which
will result in the gamma parameter being two-dimensional; the first
element is just the unconstrained criterion, the second element is log
of spread between the "unbiased" criteria. Because this model is
defined by the link function, it is also possible to fit Unequal
Variance parsimonious SDT models, or meta-d' parsimonious SDT
models. Remember that appropriate link function has to be specified
when converting between the gamma and the criteria parameters using
the gamma_to_crit function.

# NEWS (and a warning)

The bhsdtr package now handles hierarchical or non-hierarchical Equal
Variance Normal SDT models, Unequal Variance Normal SDT models, and
meta-d' models with one (not meta-d') or more criteria.

Each of the three available FOP link functions ('softmax',
'log_distance', and 'log_ratio' - see below for more details) can be
used with every model.

# WARNING

The current active branch is now called v2. The interface is slightly
changed in this branch, i.e.,

- make_stan_data and make_stan_model accept the model parameter with
possible values 'sdt' (the default), 'uvsdt', and 'metad'

- make_stan_data, make_stan_model, and gamma_to_crit functions accept
  the gamma_link parameter with possible values 'softmax' (the
  default), 'log_distance', and 'log_ratio'

- if uvsdt model is used, the (possibly hierarchical) regression
  structure for the theta parameter (= log(sd_2 / sd_1), where sd_1 =
  standard deviation of the first / noise evidence distribution = 1)
  can be specified in the same way as for the delta and gamma
  parameters

- to make the code cleaner and easier to maintain the delta, gamma,
  and theta (uvsdt model: sd_ratio = exp(theta)) parameters are now
  represented as matrices, even when they are one-dimensional, so, for
  example, when there is only one dprime / delta parameter its name is
  'delta_fixed[1,1]'.

# The importance of Flexible Order-Preserving link functions in ordinal models

Response labels such as "noise" and "signal" can be viewed as values
of a nominal scale variable, however, from the point of view of Signal
Detection Theory such variables are in fact *ordinal*. That's because
in an SDT model the response "signal" corresponds to *higher* values
of internal evidence. Moreover, once the ratings (an ordinal variable)
are introduced the problem of confounding sensitivity and bias still
exists even if we consider only one kind of responses (e.g.,
"signal"); A participant may respond "high confidence" not because the
internal evidence is high, but because, for some reason, the labels
are used differently. It is just as unrealistic to assume that the
rating scale is invariant across participants, items, or conditions as
it is to assume that the SDT decision criterion is constant. It leads
to the same kind of problem when interpreting the results - observed
differences may indicate that what is supposed to be captured by the
ratings is different, or that the way the ratings are used is
different.

Consider a typical ordinal-scale variable in psychology, such as PAS
ratings, confidence ratings, or a Likert-scale item in a
questionnaire. It is natural to assume two things about such
variables:

1. *Order invariance*, i.e., whatever latent value X this outcome is
supposed to represent, higher observed values correspond to higher
values of X, e.g., higher confidence ratings correspond to higher
value (more "signal-like") of internal evidence in an SDT model, or
higher values in a questionnaire item correspond to higher values of
the property X measured by the questionnaire. When order invariance
does not hold, it indicates that the process of generating the
responses changed in a *qualitative* way, e.g., the responses in a
binary classification task were reversed because task instructions
were misunderstood, or some of the possible responses in a
Likert-scale item were interpreted in a way that was not intended by
the authors of the questionnaire.

2. *Lack of scale invariance*, i.e., the thresholds that correspond to
the discrete outcome values or labels may differ between participants,
items, or conditions, or may covary with numerical predictors. In
fact, it would be more than just surprising if evidence was found that
the mapping between the values of ordinal responses and the latent
values captured by these responses is constant between participants,
conditions or items, since the thresholds that correspond to such
responses are parts of a psychological mechanism which is certain to
be more or less unique to each participant and cannot be assumed to be
invariant across changing conditions.

Whenever typical ordinal variables are used, there is a possibility of
confounding "response bias", which in this case corresponds to the way
the response categories are used to label e.g., some internal
experience, and the internal experience itself. This problem is seen
as important in the context of binary classification tasks and SDT
theory, but it often seems to be ignored in other contexts. For
example, in the context of IRT modelling this is known as the 'item
parameter invariance' problem and it does not seem to be studied very
frequently.

Let's use the term *Flexible Order-Preserving link function* (FOP, see
the news section for more details) to denote an isomorphic function
that maps the space of ordered real vectors (i.e., *v<sub>j</sub> >
v<sub>i</sub>* if *j > i*) to the space of unresctricted real vectors
*&gamma;* in such a way that:

1. the order is preserved in a sense that *v<sub>i</sub>* is mapped to
*&gamma;<sub>i</sub>*

2. *individual* thresholds/criteria become "free", i.e., each element
of *&gamma;* is unbounded and can be related in an arbitrary way to
nominal (e.g., participants, items, conditions) or to numerical
predictors.

*By using a FOP link function any model which represents an ordinal
variable in terms of ordered thresholds can be supplemented with a
hierarchical linear regression structure in a way that accounts for
the effects in latent values as well as for the effects in thresholds*

A model that (unrealistically) assumes that the pattern of thresholds'
placement is constant cannot account for the possibility of
response/scale bias; If all the thresholds are shifted by the same
amount in one direction the observed effects are the same as if the
thresholds stayed the same but the latent value changed. When the
thresholds can be related in a different way to various predictors
deconfounding of latent values from scale bias becomes possible. Once
we assume something about the distribution of latent values it may be
possible to estimate *non-uniform* changes in the thresholds, but only
if the model can account for such effects. FOP link functions make
many such models possible. Because ordinal models are non-linear,
supplementing them with a hierarchical linear regression structure may
solve the problem of interval and point estimate bias introduced by
aggregating the data or by otherwise ignoring hierarchical data
structure. The Bayesian hierarchical SDT model as implemented in the
bhsdtr package is only one such example.

# NEWS: Two additional link functions for the criteria + a definition of a generalized link function

In the current version of the package there are three link functions
for the SDT criteria to choose from. One is the link function
described in the preprint - this is now called "softmax". This link
function (softmax followed by inverse normal CDF) is computationally
intensive and makes the task of specifying the priors for the gamma
vector difficult.

The two new simple link functions *preserve the ordering of the
criteria and at the same time allow for individual criteria effects*,
which was arguably the main contribution of the bhsdtr package in its
previous version. The default values of the *gamma_sd* (fixed effects
specification) and the *gamma_scale* parameters (random effects
specification) for the new link functions are now set to 2, but this
is based on a small number of tests with real datasets.

Note also that adding the *init_r = l* where l < 2 argument to the
*stan* function call limits the range of initial values to (-l, l)
instead of the default range (-2, 2). This helps with all the
link functions, because a value close to 2 (or -2) is often well
outside the range of reasonable values for some of the gamma/delta
parameters.

Anyway, the unconstrained gamma vector can be mapped to the ordered
criteria vector in many useful ways. Note that the main criterion (the
K/2 threshold) considered in isolation is an unconstrained
parameter. The rest of the criteria can be represented as
log-distances between criteria or as log-ratios of distances between
criteria. For example, the K/2+1 criterion can be represented as
log(c_<sub>K+1</sub> - c<sub>K/2</sub>). This general idea leads to
some intuitive solutions. One is:

the main criterion is unconstrained:

c_<sub>K/2</sub> = &gamma;<sub>K/2</sub>

the criteria above the main criterion are represented as
log-distances, e.g.:

c<sub>K/2+3</sub> = c<sub>K/2+2</sub> + exp(&gamma;<sub>K/2+3</sub>)

and similarly for the criteria below the main criterion, e.g.:

c<sub>K/2-3</sub> = c<sub>K/2-2</sub> - exp(&gamma;<sub>K/2-3</sub>)

This is the "log_distance" gamma link function. The prior for
&gamma;<sub>K/2</sub> is now easy to specify, because this element of
the &gamma; vector represents the position of the main criterion
relative to the midpoint between the evidence distribution means,
i.e., the value of 0 corresponds to no bias and the positive
(negative) values correspond to the tendency to respond "noise"
("signal"). The priors for all the other elements of the &gamma;
vector are almost as easy to specify. For example, the assumption that
the average distance between the criteria is probably .5 can be
represented by setting the means of the priors for the &gamma; vector
(except for &gamma;<sub>K/2</sub>) at log(.5).

The other link function is called "log_ratio". The K/2th element again
represents the main criterion, the &gamma;<sub>K/2+1</sub> element
represents log(c<sub>K/2+1</sub> - c<sub>K/2</sub>), which I like to
call the "spread" parameter, because all the other distances are
represented in terms of this one. The &gamma;<sub>K/2-1</sub> element
represents the assymetry between the lower and the upper spread of the
criteria which are next to the main criterion, i.e., the following
log-ratio of distances (hence the name of the link function):
log((c<sub>K/2</sub> - c<sub>K/2-1</sub>) / (c<sub>K/2+1</sub> -
c<sub>K/2</sub>)). The elements &gamma;<sub>K/2+i</sub> where i > 1
also represent ratios of distances, i.e., &gamma;<sub>K/2+i</sub> =
log((c<sub>K/2+i</sub> - c<sub>K/2+i-1</sub>) / (c<sub>K/2+1</sub> -
c<sub>K/2</sub>)), and I like to call them "upper consistency"
parameters. The elements &gamma;<sub>K/2-i</sub> where i > 1 are
"lower consistency" parameters, i.e., &gamma;<sub>K/2-i</sub> =
log((c<sub>K/2-i+1</sub> - c<sub>K/2-i</sub>) / (c<sub>K/2</sub> -
c<sub>K/2-1</sub>)). In SDT models the reasonable prior for the
log-ratio parameters has mean = log(1) = 0.

For those who enjoy this kind of thing, here is the *generalized link
function for ordered criteria*:

1) choose an index i between 1 and K-1, this will be your
unconstrained parameter

2) represent c<sub>i</sub> as &gamma;<sub>i</sub>

3) choose an index j from the remaining K-2 indices

4) represent c<sub>j</sub> as log of distance, i.e., c<sub>j</sub> +
exp(&gamma;<sub>i</sub>) or c<sub>j</sub> - exp(&gamma;<sub>i</sub>),
depending on which threshold is supposed to be to the right of the
other

5) choose an index k from the remaining K-3 indices

6) represent c<sub>k</sub> as log of distance between c<sub>k</sub>
and c<sub>i</sub> or between c<sub>k</sub> and c<sub>j</sub> or as log
of distance between c<sub>k</sub> and c<sub>i</sub> (or c<sub>k</sub>
and c<sub>j</sub>) divided by the distance between c<sub>j</sub> and
c<sub>i</sub>

7) etc.

A broad class of models, some of which can be thought of as simplified
in a meaningful way, can be obtained just by restricting the values of
the elements of the &gamma; vector when using the two new link
functions. For example, by using the "log_ratio" link function and
fixing all the ratios at log(1) = 0 we get, as a special case, the
parsimonious SDT model as described in this great
[paper](https://link.springer.com/article/10.3758/s13428-019-01231-3)
by Selker, van den Bergh, Criss, and Wagenmakers. The &gamma; vector
can also be constrained in other ways, in particular, the constraints
can be soft (i.e., priors with small SDs) which means that a continuum
of more and more simplified models can be obtained. Moreover, the
effects of *numerical* predictors (e.g., presentation time, stimulus
intensity) on the &gamma; parameters may have a reasonably intuitive
interpretation or may not require a highly flexible polynomial to
approximate the relationship well.

In order to use the new link functions the appropriate name has to be
specified when calling the *make_stan_data*, *make_stan_model*, and
*gamma_to_crit* functions, as described in the documentation.

# bhsdtr

The bhsdtr (short for Bayesian Hierarchical Signal Detection Theory
with Ratings) package implements a novel method of Bayesian inference
for hierarchical or non-hierarchical equal variance normal Signal
Detection Theory models with one or more criteria. It uses the
state-of-the-art platform [Stan](http://mc-stan.org/) for sampling
from posterior distributions. Our method can accommodate binary
responses as well as additional ratings and an arbitrary number of
nested or crossed random grouping factors. SDT parameters can be
regressed on additional predictors within the same model via
intermediate unconstrained parameters and the model can be extended
by using automatically generated human-readable Stan code as a
template.

## Background

The equal-variance SDT with one criterion is almost (d' is not constrained to be non-negative) equivalent to probit
regression (see [this
paper](http://www.columbia.edu/~ld208/psymeth98.pdf) by DeCarlo) which
means that any software capable of fitting hierarchical generalized
linear models can be used to fit the hierarchical version of
equal-variance SDT *with one criterion and possibly negative d'*. However, the single-criterion
SDT model is untestable because the data and the model have the same
dimensionality. The main reason for using SDT is to deconfound
sensitivity and bias. This can only be achieved if an SDT model is
approximately true, but there is no way to test it in the
single-criterion case. An SDT model becomes testable (e.g., by
comparing the theoretical and the observed ROC curves) when it is
generalized - by introducing additional criteria - to the version that
accommodates ratings (e.g., "I am almost certain that this item is
new").

SDT is a *non-linear* model. An immediate consequence of non-linearity
is that inference based on data aggregated over "random" grouping
factors (such as subjects or items) is invalid because the resulting
estimates are biased (see [this
paper](http://rouder.psyc.missouri.edu/sites/default/files/morey-jmp-zROC-2008_0.pdf)
by Morey, Pratte, and Rouder for a demonstration, or see our
[preprint](http://dx.doi.org/10.23668/psycharchives.2725) for an even
more striking demonstration). The only way to avoid this problem is to
model the (possibly correlated) effects of all the relevant random
grouping factors.

A subset of hierarchical SDT models with ratings can be fitted using
hierarchical ordered regression models, such as the cumulative model
in the excellent
[brms](https://cran.r-project.org/web/packages/brms/index.html)
package. As we explain [in the
preprint](http://dx.doi.org/10.23668/psycharchives.2725), the d'
parameter is non-negative by definition and ignoring this assumption
may lead to problems if a bayesian SDT model is used. Without some
modifications (which can be done in the brms package) hierarchical
ordinal regression models do not restrict the d' to be non-negative,
because in such models d' is just the unconstrained linear regression
slope that represents the effect of the stimulus class ("noise" or
"signal"). Moreover, in typical situations, it does not make much
sense to assume that the d' random effects are normally
distributed. Finally, in the cumulative model the parameters that
correspond to the criteria in an SDT model cannot be affected
differently by the same grouping factor (i.e., the effects are
constant across categories), because the criteria in this model are
simply additional effects in the linear part. This means that the
model assumes that the pattern of the criteria is the same for every
participant (or item, etc.). Participants differ in their criteria
placement patterns and so the data from a typical rating experiment
cannot be independent given an SDT model with ratings represented as a
cumulative model.

In the bhsdtr package the generalized SDT model is supplemented with a
hierarchical linear regression structure (normally distributed
correlated random effects) thanks to a novel parametrization described
in [this preprint](http://dx.doi.org/10.23668/psycharchives.2725)
(which is now under review), and (more concisely) in the package
documentation. [Here](https://github.com/boryspaulewicz/bhsdtr/tree/master/inst/preprint/analysis_script.R)
is the annotated R script that performs all the analyses and produces
all the tables and some of the figures in the paper.

## Features

The bhsdtr package can be used to:

- fit generalized (more than one criterion), [meta-d'](http://www.columbia.edu/~bsm2105/type2sdt/), or basic (one criterion) equal variance SDT models
- fit hierarchical or non-hierarchical (e.g., single participant) models
- assess the fit using publication-ready ROC curve and combined response distribution plots with predictive intervals calculated for the chosen alpha level
- model the dependence of the SDT parameters on additional variables (e.g., task difficulty) using separate linear regression structures for the delta (d', meta-d') and gamma (criteria) parameters

### Prerequisites

A fairly up-to-date version of [R](https://www.r-project.org/) with
[the devtools
package](https://cran.r-project.org/web/packages/devtools/index.html)
already installed.

## Installing

The bhsdtr package, together will all of its dependencies, can be
installed directly from this github repository using the devtools
package:

```
devtools::install_git('git://github.com/boryspaulewicz/bhsdtr')
```

The installed package can be loaded using:

```
library(bhsdtr)
```

## Usage example

The package contains the gabor dataset


```
data(gabor)
head(gabor)
?gabor
```

To fit a hierarchical SDT model to this data we need to create some
data structures required by the stan function. This is how you can
create the combined response variable that encodes both the binary
classification decision and rating:

```
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
```

The combined responses have to be aggregated to make the sampling more
efficient. This is done using the aggregate_resoponses function, which
requires the names of the stimulus and the combined response
variables. If the data have a hierarchical structure, this structure
has to be preserved by listing the variables that cannot be collapsed
by the aggregation step. Here we list three variables that have to be
preserved: duration, id and order.

```
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
```

Finally, the fixed and random effects structure has to be specified
using lists of R model formulae. Here we assume that d' (= exp(delta))
depends on duration (a within-subject variable) and order (a
between-subject variable), but gamma (from which the criteria
parameter vector is derived) depends only on order. There is only one
random grouping factor - id - which represents the subjects. Note that
the random effects specification is a *list of lists* of model
formulate. That's because there can be more than one random grouping
factor.

```
fixed = list(delta = ~ -1 + duration:order, gamma = ~ -1 + order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
```

Now we can start sampling (note that the 'log_distance' link function
is used for the criteria instead of the detault 'softmax' and that the
*init_r* argument is added):

```
fit = stan(model_code = make_stan_model(random, gamma_link = 'log_distance'),
    data = make_stan_data(adata, fixed, random, gamma_link = 'log_distance'),
    pars = c('delta_fixed', 'gamma_fixed',
        'delta_sd_1', 'gamma_sd_1',
        'delta_random_1', 'gamma_random_1',
        'Corr_delta_1', 'Corr_gamma_1',
        ## we need counts_new for plotting
        'counts_new'),
    init_r = .5,
    iter = 8000,
    chains = 4)
```

When the make_stan_model and make_stan_data functions are called with
the optional metad=TRUE argument the meta-d' model is fitted. There
are two delta (d') parameters in the meta-d' model and so the
delta_fixed regression coefficients form a two-row matrix: the first
row represents the fixed effects for the d' parameter and the
second row represents the fixed effects for the meta-d' parameter.

Here is how you can obtain the avarege d' values per condition:

```
samples = as.data.frame(fit)
apply(exp(samples[, grep('delta_fixed', names(samples))]), 2, mean)
```

Note that exponentiation is done first and averaging second. Because
it is a non-linear transormation doing it in the other order would not
give the correct result. The average criteria for both conditions can
be obtained by calling:

```
## First condition
apply(gamma_to_crit(samples, gamma_link = 'log_distance', 1), 2, mean)
## Second condition
apply(gamma_to_crit(samples, gamma_link = 'log_distance', 2), 2, mean)
```

Note that we have to specify the correct gamma link function, since
now there are three such functions to choose from.

The model fit can be assessed using the plot_sdt_fit function, which
produces ROC curve plots ...


```
plot_sdt_fit(fit, adata, c('order', 'duration')))
```

![ROC curve](inst/preprint/roc_fit.png)

... or combined response distribution plots:

```
plot_sdt_fit(fit, adata, c('order', 'duration'), type = 'response')
```

![Combined response distributions](inst/preprint/response_fit.png)

