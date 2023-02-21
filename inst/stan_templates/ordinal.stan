// -*- coding: utf-8 -*-

data{
  int<lower=0,upper=1> PRINT;
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> Kb2;
  int<lower=0> counts[N, K];
  //cb{links$gamma=='parsimonious'}
  vector[K-1] unbiased;
  //ce
  //cb{links$gamma=='softmax'}
  real thresholds_scale;
  //ce
  //cb{model %in% c('sdt', 'uvsdt', 'metad', 'dpsdtcor', 'dpsdt')}
  vector[N] stim_sign;
  //ce
  //cb{model %in% c('dpsdtcor', 'dpsdt')}
  vector[N] modifier;
  //ce
  //cb,fpariter
  int<lower=1> PAR_size;
  int<lower=1> PAR_size_;
  //ce
  //cb,fpariter
  // Fixed effects matrices
  int<lower=1> X_PAR_ncol;
  row_vector[X_PAR_ncol] X_PAR[N];
  matrix[PAR_size, X_PAR_ncol] PAR_is_fixed;
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_value;
  //ce
  //cb,fpariter
  // Priors
  matrix[PAR_size, X_PAR_ncol] PAR_prior_fixed_mu;
  row_vector<lower=0>[X_PAR_ncol] PAR_prior_fixed_sd[PAR_size];
  //ce
  //cb{grepl('l', bounds.fe[[par]])},fpariter
  real PAR_prior_fixed_lb;
  //ce
  //cb{grepl('u', bounds.fe[[par]])},fpariter
  real PAR_prior_fixed_ub;
  //ce
  //cb,rpariter
  // Random effects matrices
  int<lower=1> PAR_group_max_G;
  int<lower=1,upper=PAR_group_max_G> PAR_group_G[N];
  int<lower=1> Z_PAR_ncol_G;
  row_vector[Z_PAR_ncol_G] Z_PAR_G[N];
  //ce
  //cb,rpariter
  // Random effects priors
  real<lower=1> PAR_prior_random_nu_G;
  matrix[PAR_size,Z_PAR_ncol_G] PAR_prior_random_scale_G;
  //ce
  //cb{grepl('l', bounds.sd[[par]][g])},rpariter
  real<lower=0> PAR_prior_random_sd_lb_G;
  //ce
  //cb{grepl('u', bounds.sd[[par]][g])},rpariter
  real<lower=0> PAR_prior_random_sd_ub_G;
  //ce
  //cb{bounds.sd[[par]][g] == ''},rpariter
  row_vector[Z_PAR_ncol_G] PAR_prior_random_scale_G[PAR_size];
  //ce
}

parameters{
  //cb{bounds.fe[[par]] == ''},fpariter
  // Fixed effects
  matrix[PAR_size, X_PAR_ncol] PAR_fixed;
  //ce
  //cb{bounds.fe[[par]] == 'l'},fpariter
  matrix<lower=PAR_prior_fixed_lb>[PAR_size, X_PAR_ncol] PAR_fixed;
  //ce
  //cb{bounds.fe[[par]] == 'u'},fpariter
  matrix<upper=PAR_prior_fixed_ub>[PAR_size, X_PAR_ncol] PAR_fixed;
  //ce
  //cb{bounds.fe[[par]] == 'lu'},fpariter
  matrix<lower=PAR_prior_fixed_lb,upper=PAR_prior_fixed_ub>[PAR_size, X_PAR_ncol] PAR_fixed;
  //ce
  //cb,rpariter
  // Random effects
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_G] L_corr_PAR_G;
  //ce
  //cb{bounds.sd[[par]][g] == 'l'},rpariter
  row_vector<lower=PAR_prior_random_sd_lb_G>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  //ce
  //cb{bounds.sd[[par]][g] == 'lu'},rpariter
  row_vector<lower=PAR_prior_random_sd_lb_G,upper=PAR_prior_random_sd_ub_G>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  //ce
  //cb,rpariter
  vector[PAR_size * Z_PAR_ncol_G] PAR_z_G[PAR_group_max_G];
  //ce
}

transformed parameters{
  //cb,fpariter
  // Fixed effects with possibly fixed values
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;
  //ce
  //cb,rpariter
  // Random effects
  matrix[PAR_size, Z_PAR_ncol_G] PAR_random_G[PAR_group_max_G];
  //ce
  //cb,rpariter
  // Vectorized matrix of random effects' standard deviations
  vector<lower=0>[PAR_size * Z_PAR_ncol_G] PAR_sd_G_v;
  //ce
  //cb,rpariter
  matrix[PAR_size * Z_PAR_ncol_G, PAR_size * Z_PAR_ncol_G] corr_PAR_G;
  //ce
  //cb,fpariter
  vector[PAR_size] PAR_fixef;
  vector[PAR_size] PAR_ranef;
  vector[PAR_size] PAR; // the group-specific effect
  vector[PAR_size_] PAR_; // = invlink(PAR)
  //ce
  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];
  //cb{model %in% c('sdt', 'uvsdt', 'dpsdtcor', 'dpsdt')}
  real shift;
  //ce
  //cb{model == 'metad'}
  real shift1;
  real shift2;
  //ce
  //cb{model == 'metad'}
  vector[2] normalization; // meta-d'
  //ce
  //cb,fpariter
  // Fixing fixed effects if requested
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      if(PAR_is_fixed[i, j] == 1){
        PAR_fixed_[i, j] = PAR_fixed_value[i, j];
      }else{
        PAR_fixed_[i, j] = PAR_fixed[i, j];
      }
  //ce
  //cb,rpariter
  // Random effects correlation matrices
  corr_PAR_G = multiply_lower_tri_self_transpose(L_corr_PAR_G); // corr_PAR_G = L_corr_PAR_G * L_corr_PAR_G';
  //ce
  //cb,rpariter
  // Vectorization of random effects' sd matrices, column major order
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G)
      PAR_sd_G_v[i + (j - 1) * PAR_size] = PAR_sd_G[i, j];
  for(g in 1:PAR_group_max_G)
    PAR_random_G[g] = to_matrix(diag_pre_multiply(PAR_sd_G_v, L_corr_PAR_G) * PAR_z_G[g], PAR_size, Z_PAR_ncol_G);
  //ce
  for(n in 1:N){
    //cb,fpariter
    PAR_fixef = PAR_fixed_ * X_PAR[n]';
    for(i in 1:PAR_size)
      PAR_ranef[i] = 0; // this is how I handle the case when there are no random effects
    //ce
    //cb,rpariter
    PAR_ranef = PAR_ranef + PAR_random_G[PAR_group_G[n]] * Z_PAR_G[n]';
    //ce
    // Applying the inverse link functions
    //cb{links$gamma == 'id_log'}
    gamma_[Kb2] = gamma_fixef[Kb2] + gamma_ranef[Kb2]; // this threshold is unconstrained
    if(K > 2){
      for(k in 1:(Kb2 - 1))
        // For instance, if K = 4 then Kb2 = 2, thr2 = gamma2 [computed above], and thr1 = thr2 - gamma_fixef[1] * exp(gamma_ranef[1])
        gamma_[Kb2 - k] = gamma_[Kb2 - k + 1] - gamma_fixef[Kb2 - k] * exp(gamma_ranef[Kb2 - k]);
      for(k in (Kb2 + 1):(K - 1))
        // ... and thr3 = thr2 + gamma_fixef[3] * exp(gamma_ranef[3])
        gamma_[k] = gamma_[k - 1] + gamma_fixef[k] * exp(gamma_ranef[k]);
    }
    //ce
    //cb{links$gamma != 'id_log'}
    gamma = gamma_fixef + gamma_ranef;    
    //ce
    //cb{links$gamma == 'identity'}
    for(k in 1:(K - 1))
      gamma_[k] = gamma[k];
    //ce
    //cb{links$gamma == 'log_distance'}
    gamma_[Kb2] = gamma[Kb2];
    if(K > 2){
      for(k in 1:(Kb2 - 1))
        gamma_[Kb2 - k] = gamma_[Kb2 - k + 1] - exp(gamma[Kb2 - k]);
      for(k in (Kb2 + 1):(K - 1))
        gamma_[k] = gamma_[k - 1] + exp(gamma[k]);
    }
    //ce
    //cb{links$gamma == 'log_ratio'}
    gamma_[Kb2] = gamma[Kb2];
    if(K > (Kb2 + 1)){
      // spread
      gamma_[Kb2 + 1] = gamma_[Kb2] + exp(gamma[Kb2 + 1]);
      if(Kb2 > 1){
        // symmetry
        gamma_[Kb2 - 1] = gamma_[Kb2] - exp(gamma[Kb2 - 1]) * (gamma_[Kb2 + 1] - gamma_[Kb2]);
        if(K > 4){
          for(k in 1:(Kb2 - 2)){
            // upper consistency
            gamma_[Kb2 + k + 1] = gamma_[Kb2 + k] + exp(gamma[Kb2 + k + 1]) * (gamma_[Kb2 + 1] - gamma_[Kb2]);
            // lower consistency
            gamma_[Kb2 - k - 1] = gamma_[Kb2 - k] - exp(gamma[Kb2 - k - 1]) * (gamma_[Kb2] - gamma_[Kb2 - 1]);
          }
        }
      }
    }
    //ce
    //cb{links$gamma == 'parsimonious'}
    for(k in 1:(K - 1)){
      gamma_[k] = gamma[1] + exp(gamma[2]) * unbiased[k];
    }
    //ce
    //cb{links$gamma == 'softmax'}
    gamma_ = thresholds_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma, 0))), K - 1));
    //ce
    //cb{links$gamma == 'twoparameter'}
    gamma_[Kb2] = gamma[1];
    if(K > 2){
      for(k in 1:(Kb2 - 1)){
        gamma_[Kb2 + k] = gamma_[Kb2 + k - 1] + exp(gamma[2]);
        gamma_[Kb2 - k] = gamma_[Kb2 - k + 1] - exp(gamma[2]);
      }
    }
    //ce
    //cb{par != 'gamma' && links[[par]] == 'identity'},fpariter
    PAR_ = PAR_fixef + PAR_ranef;
    //ce
    //cb{par != 'gamma' && links[[par]] == 'log'},fpariter
    PAR_ = exp(PAR_fixef + PAR_ranef);
    //ce
    //cb{par != 'gamma' && links[[par]] == 'logit'},fpariter
    PAR_ = inv_logit(PAR_fixef + PAR_ranef);
    //ce
    //cb{par != 'gamma' && links[[par]] %in% c('id_log', 'id_ratio_log')},fpariter
    for(i in 1:PAR_size)PAR_[i] = PAR_fixef[i] * exp(PAR_ranef[i]);
    //ce
    //cb{par != 'gamma' && links[[par]] == 'log_logit'},fpariter
    PAR_ = PAR_fixef + PAR_ranef;
    PAR_[1] = exp(PAR_[1]);
    PAR_[2] = inv_logit(PAR_[2]);
    //ce
    //cb{par == 'delta' && links[[par]] == 'log_logratio'},fpariter
    PAR_ = PAR_fixef + PAR_ranef;
    PAR_[1] = exp(PAR_[1]);
    PAR_[2] = PAR_[1] * exp(PAR_[2]);
    //ce
    //cb{par == 'delta' && links[[par]] == 'id_ratio_log'},fpariter
    PAR_ = PAR_fixef + PAR_ranef;
    PAR_[1] = PAR_[1];
    PAR_[2] = PAR_[1] * PAR_[2];
    //ce
    // Likelihood
    multinomial_cum[1] = 0;
    multinomial_cum[K + 1] = 1;
    //cb{model %in% c('sdt', 'uvsdt', 'dpsdtcor', 'dpsdt')}
    shift = -0.5 * stim_sign[n] * delta_[1];
    //ce
    //cb{model == 'metad'}
    shift1 = -0.5 * stim_sign[n] * delta_[1];
    shift2 = -0.5 * stim_sign[n] * delta_[2];
    //ce
    //cb{model %in% c('sdt', 'dpsdtcor', 'dpsdt')}
    for(k in 1:(K - 1))
      multinomial_cum[k + 1] = Phi(gamma_[k] + shift);
    //ce
    //cb{model == 'uvsdt'}
    for(k in 1:(K - 1)){
      if(stim_sign[n] > 0){
        multinomial_cum[k + 1] = Phi((gamma_[k] + shift) / theta_[1]);
      }else{
        multinomial_cum[k + 1] = Phi((gamma_[k] + shift));
      }
    }
    //ce
    //cb{model == 'metad'}
    for(k in 1:(K - 1))
      multinomial_cum[k + 1] = Phi(gamma_[k] + shift2);
    //ce
    //cb{model == 'ordinal'}
    for(k in 1:(K - 1))
      multinomial_cum[k + 1] = Phi(gamma_[k] - eta_[1]);
    //ce
    //cb{model == 'uvordinal'}
    sd_ratio = theta_[1];
    for(k in 1:(K - 1))
      multinomial_cum[k + 1] = Phi((gamma_[k] - eta_[1]) / sd_ratio);
    //ce
    //cb{model == 'metad'}
    normalization[1] = Phi(gamma_[Kb2] + shift1) / Phi(gamma_[Kb2] + shift2);
    normalization[2] = Phi(-(gamma_[Kb2] + shift1)) / Phi(-(gamma_[Kb2] + shift2));
    for(k in 1:K){
      multinomial_p[n, k] = (multinomial_cum[k + 1] - multinomial_cum[k]);
      if(k < (Kb2 + 1)){
        multinomial_p[n, k] = multinomial_p[n, k] * normalization[1];
      }else{
        multinomial_p[n, k] = multinomial_p[n, k] * normalization[2];
      }
    }
    //ce
    //cb{model != 'metad'}
    for(k in 1:K)
      multinomial_p[n, k] = multinomial_cum[k + 1] - multinomial_cum[k];
    //ce
    //cb{model == 'dpsdtcor'}
    for(k in 1:K)
      multinomial_p[n, k] = (1 - delta_[2]) * multinomial_p[n, k];
    if(stim_sign[n] > 0){
      multinomial_p[n, K] = multinomial_p[n, K] + delta_[2];
    }else{
      multinomial_p[n, 1] = multinomial_p[n, 1] + delta_[2];
    }
    //ce
    //cb{model == 'dpsdt'}
    if(modifier[n] == 1){
      for(k in 1:K)
        multinomial_p[n, k] = (1 - theta_[1]) * multinomial_p[n, k];
      if(stim_sign[n] > 0){
        multinomial_p[n, K] = multinomial_p[n, K] + theta_[1];
      }else{
        multinomial_p[n, 1] = multinomial_p[n, 1] + theta_[1];
      }
    }
    //ce
  }
}

model{
  //cb{bounds.fe[[par]] == ''},fpariter
  // Fixed effects' priors
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);
  //ce
  //cb{bounds.fe[[par]] == 'l'},fpariter
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]) -
        normal_lccdf(PAR_prior_fixed_lb | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);
  //ce
  //cb{bounds.fe[[par]] == 'u'},fpariter
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]) -
        normal_lcdf(PAR_prior_fixed_ub | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);
  //ce
  //cb{bounds.fe[[par]] == 'lu'},fpariter
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      target += normal_lpdf(PAR_fixed[i, j] | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]) -
        log_diff_exp(normal_lcdf(PAR_prior_fixed_ub | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]),
                     normal_lcdf(PAR_prior_fixed_lb | PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]));
  //ce
  //cb{!all(sdata[[sprintf('%s_prior_random_scale_%d',par,g)]] == 0)},rpariter
  // Random effects' sd priors
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G){
      target += cauchy_lpdf(PAR_sd_G[i, j] | 0, PAR_prior_random_scale_G[i, j]);
      target += -cauchy_lccdf(0 | 0, PAR_prior_random_scale_G[i, j]);
    }
  //ce
  //cb,rpariter
  target += lkj_corr_cholesky_lpdf(L_corr_PAR_G | PAR_prior_random_nu_G);
  //ce
  //cb,rpariter
  // Random effects before scaling
  for(g in 1:PAR_group_max_G)
    target += normal_lpdf(PAR_z_G[g] | 0, 1);
  //ce
  for(n in 1:N)
    target += multinomial_lpmf(counts[n] | multinomial_p[n]);
}

generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));
}
