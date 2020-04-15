## -*- coding: utf-8 -*-

## fixed = m$fixed; random = m$random; model = m$model; links = m$links

## parse.model.code = function(m){
##     bounds.fe = bounds.sd = list()
##     for(par in names(m$fixed)){
##         lb.fe.name = sprintf('%s_prior_fixed_lb', par)
##         ub.fe.name = sprintf('%s_prior_fixed_ub', par)
##         bounds.fe[[par]] = ''
##         if(!is.na(m$sdata[[lb.fe.name]]))
##             bounds.fe[[par]] = paste(bounds.fe[[par]], 'l', sep = '')
##         if(!is.na(m$sdata[[ub.fe.name]]))
##             bounds.fe[[par]] = paste(bounds.fe[[par]], 'u', sep = '')
##     }
##     for(par in names(m$random)){
##         bounds.sd[[par]] = rep('', length(m$random[[par]]))
##         for(g in 1:length(m$random[[par]])){
##             lb.sd.name = sprintf('%s_prior_random_sd_lb_%d', par, g)
##             ub.sd.name = sprintf('%s_prior_random_sd_ub_%d', par, g)
##             if(!is.na(m$sdata[[lb.sd.name]]))
##                 bounds.sd[[par]][g] = paste(bounds.sd[[par]][g], 'l', sep = '')
##             if(!is.na(m$sdata[[ub.sd.name]]))
##                 bounds.sd[[par]][g] = paste(bounds.sd[[par]][g], 'u', sep = '')
##         }
##     }
##     parsed = ''
##     lines = stan.template('ordinal_new.stan')
##     l = 1
##     while(l <= length(lines)){
##         if(grepl('//cb', lines[l])){
##             cb = l
##             for(k in (l + 1):length(lines))
##                 if(grepl('//ce', lines[k]))break
##             ce = k
##             header = lines[cb]
##             chunk = paste(lines[(cb + 1):(ce - 1)], collapse = '\n')
##             processed = process.chunk(header, chunk, m$model, m$links, m$fixed, m$random, m$sdata)
##             if(processed != "")
##                 parsed = paste(parsed, processed, sep = '\n')
##             l = k + 1
##         }else{
##             parsed = paste(parsed, lines[l], sep = '\n')
##             l = l + 1
##         }
##     }
##     parsed
## }

## ## matrix<BOUNDS> -> matrix or matrix<lower=x> etc
## parse.decl.bounds = function(chunk, lb.name, ub.name, sdata){
##     lb = sdata[[lb.name]][1]
##     ub = sdata[[ub.name]][1]
##     bounds = NULL
##     if(!is.na(lb))
##         bounds = c(bounds, sprintf('lower=%f', lb))
##     if(!is.na(ub))
##         bounds = c(bounds, sprintf('upper=%f', ub))
##     if(!is.null(bounds)){
##         bounds = sprintf('<%s>', paste(bounds, collapse = ','))
##     }else{
##         bounds = ''
##     }
##     gsub('<BOUNDS>', bounds, chunk)
## }

## ## Bounds possible only for fixed effects and random effects'
## ## std. devs.
## process.chunk = function(header, chunk, model, links, fixed, random, sdata){
##     parsed = NULL
##     test = regmatches(header, regexpr('\\{.*\\}', header))
##     if(length(test) == 0){
##         test = parse(text = 'TRUE')
##     }else{
##         test = parse(text = test)
##     }
##     process = T
##     ## if(length(test) != 0)
##     ##     process = eval(parse(text = test))
##     ## if(process){
##     if(grepl('fpariter', header)){
##         for(par in names(fixed)){
##             if(eval(test)){
##                 if(grepl('fbounds', header)){
##                     chunk.bounds = parse.decl.bounds(chunk, sprintf('%s_prior_fixed_lb', par),
##                                                      sprintf('%s_prior_fixed_ub', par), sdata)
##                 }else{
##                     chunk.bounds = chunk
##                 }
##                 parsed = c(parsed, gsub('PAR', par, chunk.bounds))
##             }
##         }
##     }else if(grepl('rpariter', header)){
##         for(par in names(random)){
##             chunk.par = gsub('PAR', par, chunk)
##             if(grepl('gpariter', header)){
##                 for(g in 1:length(random[[par]])){
##                     if(grepl('sdbounds', header)){
##                         chunk.par.bounds = parse.decl.bounds(chunk.par, sprintf('%s_prior_random_sd_lb_%d', par, g),
##                                                              sprintf('%s_prior_random_sd_ub_%d', par, g), sdata)
##                     }else{
##                         chunk.par.bounds = chunk.par
##                     }
##                     parsed = c(parsed, gsub('G', g, chunk.par.bounds))
##                 }
##             }else{
##                 parsed = c(parsed, chunk.par)
##             }
##         }
##     }else if(eval(test)){
##         parsed = chunk
##     }
##     paste(parsed, collapse = '\n')
## }

## old version
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
        code = adds(code, '  vector[K-1] unbiased;\n')
    }else if(links$gamma == 'softmax'){
        code = adds(code, '  real thresholds_scale;\n')
    }
    if(model %in% c('sdt', 'uvsdt', 'metad'))
        code = adds(code, '  vector[N] stim_sign;\n')
    for(par in pars){
        code = adds(code, gsub('PAR', par,
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
  // vector[PAR_size * Z_PAR_ncol_G] zeros_PAR_G;
  row_vector[Z_PAR_ncol_G] Z_PAR_G[N];
  // PAR random effects priors
  real<lower=1> PAR_prior_random_nu_G;
  row_vector<lower=0>[Z_PAR_ncol_G] PAR_prior_random_scale_G[PAR_size];\n')
        code = adds(code, gsub('G', g, part))
    }
    ## parameters block
    code = adds(code, '}\n\nparameters {
  // Fixed effects\n')
    for(par in names(fixed)){
        restr = ''
        if(links[[par]] == 'id_log')
            restr = '<lower=0>'
        code = adds(code, gsub('PAR', par, sprintf('  matrix%s[PAR_size, X_PAR_ncol] PAR_fixed;\n', restr)))
    }
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_G] L_corr_PAR_G;
  row_vector<lower=0>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  vector[PAR_size * Z_PAR_ncol_G] PAR_z_G[PAR_group_max_G];\n')
##  vector[PAR_size * Z_PAR_ncol_G] PAR_random_G_v[PAR_group_max_G];\n')
##  vector[PAR_size * Z_PAR_ncol_G] PAR_z_G[PAR_group_max_G];\n')
            code = adds(code, gsub('G', g, part))
        }
    ## transformed parameters block
    code = adds(code, '}\n\ntransformed parameters {\n')
    for(par in pars)
        code = adds(code, gsub('PAR', par,
'  // PAR fixed effects with possibly fixed values
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  matrix[PAR_size, Z_PAR_ncol_G] PAR_random_G[PAR_group_max_G];
  // vectorized matrix of PAR random effects\' standard deviations
  vector<lower=0>[PAR_size * Z_PAR_ncol_G] PAR_sd_G_v;
  matrix[PAR_size * Z_PAR_ncol_G, PAR_size * Z_PAR_ncol_G] corr_PAR_G;\n')
            code = adds(code, gsub('G', g, part))
        }
    for(par in pars)
        code = adds(code, gsub('PAR', par,
'  vector[PAR_size] PAR_fixef;
  vector[PAR_size] PAR;
  vector[PAR_size_] PAR_; // = invlink(PAR)\n'))
    for(par in names(random))
        code = adds(code, gsub('PAR', par, '  vector[PAR_size] PAR_ranef;\n'))
    code = adds(code,
'  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];\n')
    if(model %in% c('sdt', 'uvsdt', 'metad'))
        code = adds(code, '  real shift;\n')
    if(model == 'metad')
        code = adds(code, '  vector[2] normalization;\n')
    for(par in names(fixed))
        code = adds(code, gsub('PAR', par,
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
  // corr_PAR_G = L_corr_PAR_G * L_corr_PAR_G\';
  corr_PAR_G = multiply_lower_tri_self_transpose(L_corr_PAR_G);
  // vectorization PAR of random effects\' sd matrices, column major order
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G)
      PAR_sd_G_v[i + (j - 1) * PAR_size] = PAR_sd_G[i, j];
  for(g in 1:PAR_group_max_G)
    PAR_random_G[g] = to_matrix(diag_pre_multiply(PAR_sd_G_v, L_corr_PAR_G) * PAR_z_G[g], PAR_size, Z_PAR_ncol_G);\n')
    ## PAR_random_G[g] = to_matrix(PAR_random_G_v[g], PAR_size, Z_PAR_ncol_G);\n')
    ## PAR_random_G[g] = to_matrix(diag_pre_multiply(PAR_sd_G_v, L_corr_PAR_G) * PAR_z_G[g], PAR_size, Z_PAR_ncol_G);\n')
            code = adds(code, gsub('G', g, part))
        }
    code = adds(code, '  for(n in 1:N){\n')
    for(par in names(fixed))
        code = adds(code, gsub('PAR', par,
'    PAR_fixef = PAR_fixed_ * X_PAR[n]\';\n'))
    ## This is here for the ranef = ranef + ... part to work
    for(par in names(random))
        code = adds(code, gsub('PAR', par,
'    for(i in 1:PAR_size)
      PAR_ranef[i] = 0;\n'))
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'    PAR_ranef = PAR_ranef + PAR_random_G[PAR_group_G[n]] * Z_PAR_G[n]\';\n')
            code = adds(code, gsub('G', g, part))
        }
    ## invlink functions
    code = adds(code, '    // inverse link function\n')
    for(par in pars)
        if(par == 'gamma'){
            ## gamma parameters are link-transformed after adding the random effects, if any
            if(is.null(random[['gamma']])){
                part = '    gamma = gamma_fixef;\n'
            }else{
                part = '    gamma = gamma_fixef + gamma_ranef;\n'
            }
            code = adds(code, part)
            code = adds(code, stan.template(sprintf('gamma_link_%s.stan', links[[par]])))
        }else{
            if(links[[par]] == 'identity'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = PAR_fixef;\n'
                }else{
                    part = '    PAR_ = PAR_fixef + PAR_ranef;\n'
                }
                code = adds(code, gsub('PAR', par, part))
            }else if(links[[par]] == 'log'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = exp(PAR_fixef);\n'
                }else{
                    part = '    PAR_ = exp(PAR_fixef + PAR_ranef);\n'
                }
                code = adds(code, gsub('PAR', par, part))
            }else if(links[[par]] == 'id_log'){
                if(is.null(random[[par]])){
                    part = '    PAR_ = PAR_fixef;\n'
                }else{
                    part = '    for(i in 1:PAR_size)PAR_[i] = PAR_fixef[i] * exp(PAR_ranef[i]);\n'
                }
                code = adds(code, gsub('PAR', par, part))
            }
        }
    ## likelihood
    code = adds(code, stan.template(sprintf('likelihood_%s.stan', model)))
    ## model
    code = adds(code, '\n  }\n}\n\nmodel {\n')
    ## fixed effects' priors
    for(par in names(fixed)){
        if(links[[par]] == 'id_log'){
            prior = ' target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]) - \
normal_lccdf(0 | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);\n'
        }else{
            prior = ' target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);\n'
        }
        code = adds(code, gsub('PAR', par,
                                      adds(
'  // PAR fixed effects priors
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      ', prior)))
    }
    for(par in names(random))
        for(g in 1:length(random[[par]])){
            part = gsub('PAR', par,
'  // PAR random effects
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G){
      target += cauchy_lpdf(PAR_sd_G[i, j] | 0, PAR_prior_random_scale_G[i, j]);
      target += -cauchy_lccdf(0 | 0, PAR_prior_random_scale_G[i, j]);
    }
  target += lkj_corr_cholesky_lpdf(L_corr_PAR_G | PAR_prior_random_nu_G);
  for(g in 1:PAR_group_max_G)
    target += normal_lpdf(PAR_z_G[g] | 0, 1);\n')
##    target += multi_normal_cholesky_lpdf(PAR_random_G_v[g] | zeros_PAR_G, diag_pre_multiply(PAR_sd_G_v, L_corr_PAR_G));\n')
##    target += normal_lpdf(PAR_z_G[g] | 0, 1);\n')
            code = adds(code, gsub('G', g, part))
        }
if(only_prior){
    code = adds(code, '\n}\n')
}else{
   code = adds(code,
'  for(n in 1:N)
    target += multinomial_lpmf(counts[n] | multinomial_p[n]);\n}\n\n')
}
    ## generated quantities
    code = adds(code,
'generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));\n}\n')
    code
}

adds = function(...){
    strs = list(...)
    res = ''
    for(s in strs)
        res = paste(res, s, sep = '')
    res
}

stan.template = function(fname){
    paste(readLines(stan.file(fname)), collapse = '\n')
}

stan.file = function(fname){
    sprintf('%s/stan_templates/%s', path.package('bhsdtr2'), fname)
}
