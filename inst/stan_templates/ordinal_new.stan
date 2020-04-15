// -*- coding: utf-8 -*-

data{
  int<lower=0,upper=1> PRINT;
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> Kb2;
  int<lower=0> counts[N, K];
  //cb{model=='parsimonious'}
  vector[K-1] unbiased;
  //ce
  //cb{links$gamma=='softmax'}
  real thresholds_scale;
  //ce
  //cb{model %in% c('sdt', 'uvsdt', 'metad')}
  vector[N] stim_sign;
  //ce
  //cb,fpariter
  int<lower=1> PAR_size;
  int<lower=1> PAR_size_;
  //ce
  // Fixed effects matrices
  //cb,fpariter
  int<lower=1> X_PAR_ncol;
  row_vector[X_PAR_ncol] X_PAR[N];
  matrix[PAR_size, X_PAR_ncol] PAR_is_fixed;
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_value;
  //ce
  // Priors
  //cb,fpariter
  matrix[PAR_size, X_PAR_ncol] PAR_prior_fixed_mu;
  row_vector<lower=0>[X_PAR_ncol] PAR_prior_fixed_sd[PAR_size];
  //ce
  //cb{grepl('l', bounds.fe[[par]])},fpariter
  real PAR_prior_fixed_lb;
  //ce
  //cb{grepl('u', bounds.fe[[par]])},fpariter
  real PAR_prior_fixed_ub;
  //ce
  // Random effects matrices
  //cb,rpariter,gpariter
  int<lower=1> PAR_group_max_G;
  int<lower=1,upper=PAR_group_max_G> PAR_group_G[N];
  int<lower=1> Z_PAR_ncol_G;
  row_vector[Z_PAR_ncol_G] Z_PAR_G[N];
  //ce
  // Random effects priors
  //cb,rpariter,gpariter
  real<lower=1> PAR_prior_random_nu_G;
  //ce
  //cb{grepl('l', bounds.sd[[par]][g])},rpariter,gpariter
  real<lower=0> PAR_prior_random_sd_lb_G;
  //ce
  //cb{grepl('u', bounds.sd[[par]][g])},rpariter,gpariter
  real<lower=0> PAR_prior_random_sd_ub_G;
  //ce
  //cb{bounds.sd[[par]][g] == ''},rpariter,gpariter
  row_vector[Z_PAR_ncol_G] PAR_prior_random_scale_G[PAR_size];
  //ce
}

parameters{
  // Fixed effects
  //cb{bounds.fe[[par]] == ''},fpariter
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
  // Random effects
  //cb,rpariter,gpariter
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_G] L_corr_PAR_G;
  //ce
  //cb{bounds.sd[[par]][g] == 'l'},rpariter,gpariter
  row_vector<lower=PAR_prior_random_sd_lb_G>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  //ce
  //cb{bounds.sd[[par]][g] == 'lu'},rpariter,gpariter
  row_vector<lower=PAR_prior_random_sd_lb_G,upper=PAR_prior_random_sd_ub_G>[Z_PAR_ncol_G] PAR_sd_G[PAR_size];
  //ce
  //cb,rpariter,gpariter
  vector[PAR_size * Z_PAR_ncol_G] PAR_z_G[PAR_group_max_G];
  //ce
}

transformed parameters{
  // Fixed effects with possibly fixed values
  //cb,fpariter
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;
  //ce
  // Random effects
  //cb,rpariter,gpariter
  matrix[PAR_size, Z_PAR_ncol_G] PAR_random_G[PAR_group_max_G];
  //ce
  // vectorized matrix of random effects' standard deviations
  //cb,rpariter,gpariter
  vector<lower=0>[PAR_size * Z_PAR_ncol_G] PAR_sd_G_v;
  //ce
  //cb,rpariter,gpariter
  matrix[PAR_size * Z_PAR_ncol_G, PAR_size * Z_PAR_ncol_G] corr_PAR_G;
  //ce
  //cb,fpariter
  vector[PAR_size] PAR_fixef;
  vector[PAR_size] PAR_ranef;
  vector[PAR_size] PAR;
  vector[PAR_size_] PAR_; // = invlink(PAR)
  //ce
  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];
  //cb{model %in% c('sdt', 'uvsdt')}
  real shift;
  //ce
  //cb{model == 'metad'}
  real shift1;
  real shift2;
  //ce
  //cb{model == 'metad'}
  vector[2] normalization; // meta-d'
  //ce
  // Fixing fixed effects if requested
  //cb,fpariter
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      if(PAR_is_fixed[i, j] == 1){
        PAR_fixed_[i, j] = PAR_fixed_value[i, j];
      }else{
        PAR_fixed_[i, j] = PAR_fixed[i, j];
      }
  //ce
  // Random effects correlation matrices
  //cb,rpariter,gpariter
  corr_PAR_G = multiply_lower_tri_self_transpose(L_corr_PAR_G); // corr_PAR_G = L_corr_PAR_G * L_corr_PAR_G';
  //ce
  // Vectorization of random effects' sd matrices, column major order
  //cb,rpariter,gpariter
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
      PAR_ranef[i] = 0;
    //ce
    //cb,rpariter,gpariter
    PAR_ranef = PAR_ranef + PAR_random_G[PAR_group_G[n]] * Z_PAR_G[n]';
    //ce
    // Applying the inverse link functions
    gamma = gamma_fixef;
    gamma = gamma_fixef + gamma_ranef;    
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
    //cb{par != 'gamma' && links[[par]] == 'id_log'},fpariter
    for(i in 1:PAR_size)PAR_[i] = PAR_fixef[i] * exp(PAR_ranef[i]);
    //ce
    // Likelihood
    multinomial_cum[1] = 0;
    multinomial_cum[K + 1] = 1;
    //cb{model %in% c('sdt', 'uvsdt')}
    shift = -0.5 * stim_sign[n] * delta_[1];
    //ce
    //cb{model == 'metad'}
    shift1 = -0.5 * stim_sign[n] * delta_[1];
    shift2 = -0.5 * stim_sign[n] * delta_[2];
    //ce
    //cb{model == 'sdt'}
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
  }
}

model{
  // Fixed effects' priors
  //cb{bounds.fe[[par]] == ''},fpariter
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
  // Random effects' priors
  //cb,rpariter,gpariter
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_G){
      target += cauchy_lpdf(PAR_sd_G[i, j] | 0, PAR_prior_random_scale_G[i, j]);
      target += -cauchy_lccdf(0 | 0, PAR_prior_random_scale_G[i, j]);
    }
  target += lkj_corr_cholesky_lpdf(L_corr_PAR_G | PAR_prior_random_nu_G);
  //ce
  // Random effects before scaling
  //cb,rpariter,gpariter
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
