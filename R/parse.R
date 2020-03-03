## -*- coding: utf-8 -*-

## fixed = m$fixed; random = m$random; model = m$model; links = m$links

make.model.code = function(model, fixed, random, links, only_prior = F){
    pars = unique(names(fixed), names(random))
    ## data block
    code =
'data {
  int<lower=0,upper=1> PRINT;
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> Kb2;
  int<lower=0> counts[N, K];\n'
    if(links$gamma == 'parsimonious'){
        code = add.strings(code, '  vector[K-1] unbiased;\n')
    }else if(links$gamma == 'softmax'){
        code = add.strings(code, '  real thresholds_scale;\n')
    }
    if(model %in% c('sdt', 'uvsdt', 'metad'))
        code = add.strings(code, '  vector[N] stim_sign;\n')
    for(par in pars){
        code = add.strings(code, gsub('PAR', par,
'  int<lower=1> PAR_size;
  int<lower=1> PAR_size_;
  // Fixed effects matrices
  int<lower=1> X_PAR_ncol;
  row_vector[X_PAR_ncol] X_PAR[N];
  matrix[PAR_size, X_PAR_ncol] PAR_is_fixed;
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_value;
  // Priors
  matrix[PAR_size, X_PAR_ncol] PAR_prior_fixed_mu;
  row_vector<lower=0>[X_PAR_ncol] PAR_prior_fixed_sd[PAR_size];\n'
  ))}
    for(par in names(random)){
        for(g in 1:length(random[[par]]))
            part = gsub('PAR', par,
'  // PAR random effects matrices
  int<lower=1> PAR_group_max_G;
  int<lower=1,upper=PAR_group_max_G> PAR_group_G[N];
  int<lower=1> Z_PAR_ncol_G;
  row_vector[Z_PAR_ncol_G] Z_PAR_G[N];
  // PAR random effects priors
  real<lower=1> PAR_prior_nu_G;
  row_vector<lower=0>[Z_PAR_ncol_G] PAR_prior_scale_G[PAR_size];\n')
        code = add.strings(code, gsub('G', g, part))
    }
    code = add.strings(code, '}\n\n')
    ## parameters block
    code = add.strings(code, 'parameters {
  // Fixed effects\n')
    for(par in names(fixed))
        code = add.strings(code, gsub('PAR', par, '  matrix[PAR_size, X_PAR_ncol] PAR_fixed;\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_G] L_corr_PAR_G;
  row_vector<lower=0>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  vector[PAR_size * Z_PAR_ncol_G] PAR_z_G[PAR_group_max_G];\n')
            code = add.strings(code, gsub('G', g, part))
        }
    code = add.strings(code, '}\n\n')
    ## transformed parameters block
    code = add.strings(code, 'transformed parameters {\n')
    for(par in pars)
        code = add.strings(code, gsub('PAR', par,
'  // PAR fixed effects with possibly fixed values
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  matrix[PAR_size, Z_PAR_ncol_G] PAR_random_G[PAR_group_max_G];
  // vectorized matrix of PAR random effects\' standard deviations
  vector<lower=0>[PAR_size * Z_PAR_ncol_G] PAR_sd_G_;
  matrix[PAR_size * Z_PAR_ncol_G, PAR_size * Z_PAR_ncol_G] corr_PAR_G;\n')
            code = add.strings(code, gsub('G', g, part))
        }
    for(par in pars)
        code = add.strings(code, gsub('PAR', par,
'  vector[PAR_size] PAR_fixef;
  vector[PAR_size] PAR;
  vector[PAR_size_] PAR_; // = invlink(PAR)\n'))
    for(par in names(random))
        code = add.strings(code, gsub('PAR', par, '  vector[PAR_size] PAR_ranef;\n'))
    code = add.strings(code,
'  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];\n')
    if(model %in% c('sdt', 'uvsdt', 'metad'))
        code = add.strings(code, '  real shift;\n')
    if(model == 'metad')
        code = add.strings(code, '  vector[2] normalization;\n')
    for(par in names(fixed))
        code = add.strings(code, gsub('PAR', par,
'  // Fixing PAR fixed effects if requested
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      if(PAR_is_fixed[i, j] == 1){
        PAR_fixed_[i, j] = PAR_fixed_value[i, j];
      }else{
        PAR_fixed_[i, j] = PAR_fixed[i, j];
      }\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  corr_PAR_G = L_corr_PAR_G * L_corr_PAR_G\';
  // vectorization PAR of random effects\' sd matrices, column major order
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G)
      PAR_sd_G_[i + (j - 1) * PAR_size] = PAR_sd_G[i, j];
  for(g in 1:PAR_group_max_G)
    PAR_random_G[g] = to_matrix(diag_pre_multiply(PAR_sd_G_, L_corr_PAR_G) * PAR_z_G[g], PAR_size, Z_PAR_ncol_G);\n')
            code = add.strings(code, gsub('G', g, part))
        }
    code = add.strings(code, '  for(n in 1:N){\n')
    for(par in names(fixed))
        code = add.strings(code, gsub('PAR', par,
'    PAR_fixef = PAR_fixed_ * X_PAR[n]\';\n'))
    ## This is here for the ranef = ranef + ... part to work
    for(par in names(random))
        code = add.strings(code, gsub('PAR', par,
'    for(i in 1:PAR_size)
      PAR_ranef[i] = 0;\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'    PAR_ranef = PAR_ranef + PAR_random_G[PAR_group_G[n]] * Z_PAR_G[n]\';\n')
            code = add.strings(code, gsub('G', g, part))
        }
    ## invlink functions
    code = add.strings(code, '    // inverse link function\n')
    for(par in pars)
        if(par == 'gamma'){
            if(is.null(random[['gamma']])){
                part = '    gamma = gamma_fixef;\n'
            }else{
                part = '    gamma = gamma_fixef + gamma_ranef;\n'
            }
            code = add.strings(code, part)
            code = add.strings(code, stan.template(sprintf('gamma_link_%s.stan', links[[par]])))
        }else{
            if(links[[par]] == 'identity'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = PAR_fixef;\n'
                }else{
                    part = '    PAR_ = PAR_fixef + PAR_ranef;\n'
                }
                code = add.strings(code, gsub('PAR', par, part))
            }else if(links[[par]] == 'log'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = exp(PAR_fixef);\n'
                }else{
                    part = '    PAR_ = exp(PAR_fixef + PAR_ranef);\n'
                }
                code = add.strings(code, gsub('PAR', par, part))
            }else if(links[[par]] == 'id_log'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = PAR_fixef;\n'
                }else{
                    part = '    for(i in 1:PAR_size)PAR_[i] = PAR_fixef[i] * exp(PAR_ranef[i]);\n'
                }
                code = add.strings(code, gsub('PAR', par, part))
            }
        }
    ## likelihood
    code = add.strings(code, stan.template(sprintf('likelihood_%s.stan', model)))
    code = add.strings(code, '\n  }\n}\n\n')
    ## model
    code = add.strings(code, 'model {\n')
    for(par in names(fixed)){
        if(links[[par]] == 'id_log'){
            prior = ' ~ normal(PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j])T[0,];\n'
        }else{
            prior = ' ~ normal(PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);\n'
        }
        code = add.strings(code, gsub('PAR', par,
                                      add.strings(
'  // PAR fixed effects priors
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      PAR_fixed[i, j]', prior)))
    }
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G)
      PAR_sd_G[i, j] ~ cauchy(0, PAR_prior_scale_G[i, j]);
  L_corr_PAR_G ~ lkj_corr_cholesky(PAR_prior_nu_G);
  for(g in 1:PAR_group_max_G)
    PAR_z_G[g] ~ normal(0, 1);\n')
            code = add.strings(code, gsub('G', g, part))
        }
if(only_prior){
    code = add.strings(code, '\n}\n')
}else{
   code = add.strings(code,
'  for(n in 1:N)
    counts[n] ~ multinomial(multinomial_p[n]);\n}\n\n')
}
    ## generated quantities
    code = add.strings(code,
'generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));\n}\n')
    code
}

add.strings = function(string1, string2, collapse = '')
    paste(c(string1, string2), collapse = collapse)

stan.template = function(fname){
    paste(readLines(stan.file(fname)), collapse = '\n')
}

stan.file = function(fname){
    sprintf('%s/stan_templates/%s', path.package('bhsdtr2'), fname)
}

######################################################################
## not used anymore

## Parsing functions (model code from template)
parse.PAR = function(lines, par_types){
    parsed = NULL
    for(l in lines){
        if(rmatch('PAR', l)){
            for(par_type in par_types)
                parsed[length(parsed) + 1] = gsub('PAR', par_type, l)
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.likelihood = function(lines, model){
    parsed = NULL
    for(l in lines){
        if(rmatch('//likelihood', l)){
            parsed = c(parsed, readLines(stan.file(sprintf('likelihood_%s.stan', model))))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.link = function(lines, par, link){
    parsed = NULL
    for(l in lines){
        if(rmatch(sprintf('//%s-link', par), l)){
            parsed = c(parsed, readLines(stan.file(sprintf('%s_link_%s.stan', par, link))))
        }else{
            parsed = c(parsed, l)
        }
    }
    parsed
}

parse.random = function(lines, random, par_types){
    parsed = NULL
    for(part in lines){## If this is part of the random effects' specification ...
        if(rmatch(sprintf('//random-(%s)', paste(par_types, collapse = '|')), part)){
            ## ... and there are random effects in the model ...
            if(length(random) > 0)
                for(par in names(random)){
                    for(g in 1:length(random[[par]])){
                        if(!is.null(random[[par]][[g]]) & rmatch(sprintf('//random-%s', par), part))
                            parsed[length(parsed) + 1] = gsub('%', g, part)
                    }
                }
        }else{
            ## If this is not a part of random effects structure then do not change it
            parsed = c(parsed, part)
        }
        ## if there are no random effects in the model and the
        ## template part defines part of the random effects structure then ignore this part
    }
    parsed
}
