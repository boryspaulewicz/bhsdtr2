# bhsdtr2

The bhsdtr2 (short for Bayesian Hierarchical Signal Detection Theory
with Ratings version 2) package implements a novel method of Bayesian
inference for hierarchical or non-hierarchical ordinal models with
ordered thresholds, such as Equal-Variance Normal SDT,
Unequal-Variance Normal SDT, Equal Variance meta-d', Unequal Variance
meta-d', or the parsimonious (Equal- or Unequal-Variance) SDT and
meta-d' models (see this
[paper](https://link.springer.com/article/10.3758/s13428-019-01231-3)
by Selker, van den Bergh, Criss, and Wagenmakers for an explanation of
the term 'parsimonious' in this context).

The bhsdtr2 package uses the state-of-the-art platform
[Stan](http://mc-stan.org/) for sampling from posterior
distributions. The models can accommodate binary responses as well as
ordered polytomous responses and an arbitrary number of nested or
crossed random grouping factors. The parameters (e.g., d', criterion,
thresholds, ratio of standard deviations, latent mean) can be
regressed on additional predictors within the same model via
intermediate unconstrained parameters. The models can be extended by
modifying automatically generated human-readable Stan code.

## Background

Ordinal models are *non-linear*. An immediate consequence of
non-linearity is that inference based on data aggregated over grouping
factors (such as participants or items) is invalid because the
resulting estimates of all the model parameters are asymptotically
biased, often severely so (see [this
paper](http://rouder.psyc.missouri.edu/sites/default/files/morey-jmp-zROC-2008_0.pdf)
by Morey, Pratte, and Rouder for a demonstration, or see our
[preprint](http://dx.doi.org/10.23668/psycharchives.2725) for an even
more striking demonstration). The only correct solution to this
problem is to model the (possibly correlated) "random" effects of all
the relevant grouping factors.

In the bhsdtr2 package ordinal models are supplemented with a
hierarchical linear regression structure (normally distributed
correlated random effects) thanks to a novel parametrization described
in [this preprint](http://dx.doi.org/10.23668/psycharchives.2725)
which was recently accepted for publication in Behavior Research
Methods, and - more concisely - in the package documentation.

The main advantage of the bhsdtr2 (and bhsdtr) package over other
available methods of fitting ordinal models with ordered thresholds
has to do with the order-preserving link functions which are used for
the thresholds (criteria). To my knowledge, at present bhsdtr and
bhsdtr2 are the only correct available implementations of hierarchical
ordinal models with ordered thresholds, including hierarchical
SDT-like models, because both packages allow for variability in d' (or
latent mean) and in individual thresholds while respeting the
assumptions of non-negativity (the d' parameter in SDT-like models)
and order (thresholds) (see the
[preprint](http://dx.doi.org/10.23668/psycharchives.2725) for more
details).

Without the order-preserving link functions it is impossible to
correctly model the effects in *individual* thresholds, including the
possibly ubiquitous individual differences (i.e., participant effects)
in the *pattern* of threshold placement. Note that the preprint covers
only the SDT models and only one order-preserving link function,
whereas there are now five such functions (softmax, log_distance,
log_ratio, twoparamter and parsimonious) to choose from in bhsdtr and
bhsdtr2.

### Prerequisites

A fairly up-to-date version of [R](https://www.r-project.org/) with
[the devtools
package](https://cran.r-project.org/web/packages/devtools/index.html)
already installed.

## Installing the package

The bhsdtr2 package, together will all of its dependencies, can be
installed directly from this github repository using the devtools
package:

```
devtools::install_git('git://github.com/boryspaulewicz/bhsdtr2')
```

## Usage example

The package contains the gabor dataset

```
library(bhsdtr2)
library(rstan)
data(gabor)
head(gabor)
?gabor
```

To fit a hierarchical SDT model with multiple thresholds to this data
we need to create the combined response variable that encodes both the
binary classification decision and rating:

```
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
```

The combined responses have to be aggregated to make the sampling more
efficient. This is done internally by the bhsdtr function. Here is how
you can fit the hierarchical Equal Variance Normal SDT model, in which
we assume that d' depends on duration (a within-subject variable) and
order (a between-subject variable), the effect of duration may vary
between the participants, and the thresholds, which may also vary
between the particupants, depend only on order:

```
(m = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id)),
            r ~ stim,
            gabor))
```

Here is how you can fit the Unequal Variance version of this model, in
which the ratio of the standard deviations may vary between the
participants:

```
(m = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id), sdratio ~ 1 + (1 | id)),
            r ~ stim,
            gabor))
```

Here is how you can fit the hierarchical meta-d' model (just replace
'dprim' with 'metad' in the first model formula):

```
(m = bhsdtr(c(metad ~ duration * order + (duration | id), thr ~ order + (1 | id))
            r ~ stim,
            gabor))
```

or the hierarchical Unequal Variance meta-d' model:

```
(m = bhsdtr(c(metad ~ duration * order + (duration | id), thr ~ order + (1 | id), sdratio ~ 1 + (1 | id)),
            r ~ stim,
            gabor))
```

etc.

By default, the bhsdtr function uses the log link function for the d'
parameter and the log_distance link function for the thresholds. If
not asked to do otherwise, it quickly fits the model using maximum
likelihood optimization. To use the stan sampler the fit_method =
'stan' argument has to be added to the function call.

```
(m = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id)),
            r ~ stim,
            gabor, fit_method = 'stan'))
```

Here is how you can obtain the point estimates (ML fit) or posterior
samples of d' (or any other parameter) for every condition:

```
samples(m, 'dprim')
```

Even though the model was parametrized in such a way that the fixed
effects of duration, order and their interaction represented
differences in delta (the unconstrained version of the d' parameter),
the samples function correctly recovers posterior d' samples for every
condition by doing some tricks with model matrices, so that you don't
have to.

Here is how you can obtain the group-specific (here
participant-specific) posterior samples for every condition:

```
samples(m, 'dprim', group = 1)
```

The group = 1 argument indicates that we are interested in the effects
represented by the first (here the only) random effects model
formula (1 | id).

The model fit can be assessed using the plot_sdt_fit function, which
produces ROC curve plots ...


```
plot(m, adata, c('order', 'duration'), type = 'roc')
```

![ROC curve](inst/plots/roc_fit.png)

... or (conditional) response distribution plots (the default):

```
plot(m, c('order', 'duration'))
```

![Combined response distributions](inst/plot/response_fit.png)

# The importance of order-preserving link functions in ordinal models

Response labels such as "noise" and "signal" can be viewed as values
of a nominal scale variable, however, from the point of view of Signal
Detection Theory such variables are in fact *ordinal*. That's because
in an SDT model the response "signal" corresponds to *higher* values
of internal evidence. Moreover, once the ratings (an ordinal variable)
are introduced the problem of confounding sensitivity and bias still
exists even if we consider only one kind of responses (e.g., only
"signal" responses); A participant may respond "high confidence" not
because the internal evidence is high, but because, for some reason,
the labels are used differently. It is just as unrealistic to assume
that the rating scale is invariant across participants, items, or
conditions as it is to assume that the SDT decision criterion is
constant. It leads to the same kind of problem when interpreting the
results - observed differences may indicate that what is supposed to
be captured by the ratings is different, or that the way the ratings
are used (the average position and pattern of the thresholds) is
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
the response categories may differ between participants, items, or
conditions, or may covary with quantitative predictors. In fact, it
would be more than just surprising if evidence was found that the
mapping between the values of ordinal responses and the latent values
captured by these responses is constant between participants,
conditions or items, since the thresholds that correspond to such
responses are parts of the psychological mechanism which is certain to
be more or less unique to each participant and cannot be assumed to be
invariant across changing conditions.

Whenever typical ordinal variables are used, there is a possibility of
confounding "response bias", which in this case corresponds to the way
the response categories are used to label e.g., some internal
experience, and the internal experience itself. This problem is seen
as important in the context of binary classification tasks and SDT
theory, but it often seems to be ignored in other contexts. For
example, in the context of IRT modelling this is known as the 'item
parameter invariance' problem and when it is studied it is usually
reduced to Differential Item Functioning (DIF). However, DIF is a
population-level effect and in no way captures the individual
differences in the *pattern* of the thresholds.

An *order-preserving link function* is an isomorphic function that
maps the space of ordered real vectors (i.e., *v<sub>j</sub> >
v<sub>i</sub>* if *j > i*) to the space of unresctricted real vectors
*&gamma;* in such a way that:

1. the order is preserved in a sense that *v<sub>i</sub>* is mapped to
*&gamma;<sub>i</sub>*

2. *individual* elements (e.g., thresholds) become "free", i.e., each
element of *&gamma;* is unbounded and can be related in an arbitrary
way to nominal (e.g., participants, items, conditions) or to
quantitative predictors.

*By using an order-preserving link function any model which represents
an ordinal variable in terms of ordered thresholds can be supplemented
with a hierarchical linear regression structure in a way that accounts
for the effects in the latent values as well as for the effects in the
thresholds*

A model that (usually unrealistically) assumes that the pattern of
thresholds' placement is constant across participants or conditions
cannot account for the possibility of response/scale bias; If all the
thresholds are shifted by the same amount in one direction the
observed effects are the same as if the thresholds stayed the same but
the latent value changed. It is only when the thresholds and the
latent values are assumed to be related to different variables
(selective influence) that deconfounding of latent values from scale
bias becomes possible. Order-preserving link functions make many such
models possible. Because ordinal models are non-linear, supplementing
them with a hierarchical linear regression structure may solve the
problem of asymptotic interval and point estimate bias introduced by
aggregating the data or by otherwise ignoring hierarchical data
structure.

# Description of order-preserving link functions

In the current version of the package there are five link functions
for the thresholds to choose from. One is the link function described
in the preprint - this is now called "softmax". This link function
(softmax followed by inverse normal CDF) is quite complicated and
makes the task of specifying the priors for the gamma vector
difficult.

The two new simple link functions also preserve the ordering of the
thresholds and at the same time allow for individual threshold
effects, which was arguably the main contribution of the bhsdtr
package in its previous version.

The unconstrained gamma vector can be mapped to the ordered thresholds
vector in many useful ways. Note that the middle threshold (the K/2th
threshold) considered in isolation is an unconstrained parameter. The
rest of the thresholds can be represented as log-distances between
thresholds or as log-ratios of distances between thresholds. For
example, the K/2+1th threshold can be represented as
log(c_<sub>K+1</sub> - c<sub>K/2</sub>). This general idea leads to
some intuitive solutions. One is:

the middle threshold is unconstrained:

c_<sub>K/2</sub> = &gamma;<sub>K/2</sub>

the thresholds above the main threshold are represented as
log-distances, e.g.:

c<sub>K/2+3</sub> = c<sub>K/2+2</sub> + exp(&gamma;<sub>K/2+3</sub>)

and similarly for the thresholds below the main threshold, e.g.:

c<sub>K/2-3</sub> = c<sub>K/2-2</sub> - exp(&gamma;<sub>K/2-3</sub>)

This is the "log_distance" gamma/threshold link function. The prior
for &gamma;<sub>K/2</sub> is now easy to specify, because this element
of the &gamma; vector represents the position of the main threshold
relative to the midpoint between the evidence distribution means,
i.e., the value of 0 corresponds to no bias and the positive
(negative) values correspond to the tendency to respond "noise"
("signal"). The priors for all the other elements of the &gamma;
vector are almost as easy to specify. For example, the assumption that
the average distance between the thresholds is probably .5 can be
represented by setting the means of the priors for the &gamma; vector
(except for &gamma;<sub>K/2</sub>) at log(.5).

The other link function is called "log_ratio". The K/2th element again
represents the main threshold, the &gamma;<sub>K/2+1</sub> element
represents log(c<sub>K/2+1</sub> - c<sub>K/2</sub>), which I like to
call the "spread" parameter, because all the other distances are
represented in terms of this one. The &gamma;<sub>K/2-1</sub> element
represents the assymetry between the lower and the upper spread of the
thresholds which are next to the main threshold, i.e., the following
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

For those who enjoy this kind of thing, here is the *generalized* link
function for ordered thresholds:

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
