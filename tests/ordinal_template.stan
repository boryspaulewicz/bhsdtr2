// -*- coding: utf-8 -*-

// Every kind of parameter is represented using the type vector, and
// so the effects are represented by a matrix, which may have more
// than one row (if the vector of parameters has length > 1).

data {
  int<lower=0,upper=1> PRINT;
  int<lower=1> N;
  int<lower=2> K;
  // Kb2 = K/2 is here to avoid the (irrelevant) warning about integer
  // division
  int<lower=1> Kb2;
  // for the parsimonious link function
  vector[K-1] unbiased;
  real thresholds_scale;
  // used only in sdt-like models
  vector[N] stim_sign;
  int<lower=0> counts[N, K];
  // this will be replaced with delta_size etc.
  int<lower=1> PAR_size;
  int<lower=1> PAR_size_;
  int<lower=1> X_PAR_ncol;
  row_vector[X_PAR_ncol] X_PAR[N];
  matrix[PAR_size, X_PAR_ncol] PAR_is_fixed;
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_value;
  int<lower=1> PAR_group_max_%; //random-PAR
  int<lower=1,upper=PAR_group_max_%> PAR_group_%[N]; //random-PAR
  int<lower=1> Z_PAR_ncol_%; //random-PAR
  row_vector[Z_PAR_ncol_%] Z_PAR_%[N]; //random-PAR
  // Priors
  matrix[PAR_size, X_PAR_ncol] PAR_prior_fixed_mu;
  row_vector<lower=0>[X_PAR_ncol] PAR_prior_fixed_sd[PAR_size];
  real<lower=1> PAR_prior_nu_%; //random-PAR
  row_vector<lower=0>[Z_PAR_ncol_%] PAR_prior_scale_%[PAR_size]; //random-PAR
}

parameters {
  matrix[PAR_size, X_PAR_ncol] PAR_fixed;
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_%] L_corr_PAR_%; //random-PAR
  row_vector<lower=0>[Z_PAR_ncol_%] PAR_sd_%[PAR_size]; //random-PAR
  vector[PAR_size * Z_PAR_ncol_%] PAR_z_%[PAR_group_max_%]; //random-PAR
}

transformed parameters {
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;
  matrix[PAR_size, Z_PAR_ncol_%] PAR_random_%[PAR_group_max_%]; //random-PAR
  // vectorized matrix of random effects' standard deviations
  vector<lower=0>[PAR_size * Z_PAR_ncol_%] PAR_sd_%_; //random-PAR
  matrix[PAR_size * Z_PAR_ncol_%, PAR_size * Z_PAR_ncol_%] corr_PAR_%; //random-PAR
  vector[PAR_size] PAR;
  vector[PAR_size_] PAR_; // = invlink(par)
  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];
  // used only in the sdt model family
  real shift;
  // used only in the metad model
  vector[2] normalization;

  // fixing fixed effects if requested
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      if(PAR_is_fixed[i, j] == 1){
        PAR_fixed_[i, j] = PAR_fixed_value[i, j];
      }else{
        PAR_fixed_[i, j] = PAR_fixed[i, j];
      }

  // PAR random effects
  corr_PAR_% = L_corr_PAR_% * L_corr_PAR_%';
  // vectorization PAR of random effects' sd matrices, column major order
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_%)
      PAR_sd_%_[i + (j - 1) * PAR_size] = PAR_sd_%[i, j];
  for(g in 1:PAR_group_max_%)
    PAR_random_%[g] = to_matrix(diag_pre_multiply(PAR_sd_%_, L_corr_PAR_%) * PAR_z_%[g], PAR_size, Z_PAR_ncol_%);

  if(PRINT == 1){
    print("PRIORS: ");
    print("PAR_prior_fixed_mu "); for(i in 1:PAR_size)print(PAR_prior_fixed_mu[i,]);
    print("PAR_prior_fixed_sd"); for(i in 1:PAR_size)print(PAR_prior_fixed_sd[i,]);
    print("PAR_prior_scale_%"); for(i in 1:PAR_size)print(PAR_prior_scale_%[i,]); //random-PAR
    print("PAR_prior_nu_%"); print(PAR_prior_nu_%); //random-PAR
    print("INITIAL VALUES: ");
    print("PAR_fixed"); for(i in 1:PAR_size)print(PAR_fixed[i,]);
    for(g in 1:PAR_group_max_%)print("PAR_z_%[", g, "] = ", PAR_z_%[g]); //random-PAR
    print("PAR_sd_% = ", PAR_sd_%); //random-PAR
  }

  for(n in 1:N){
    PAR_fixef = PAR_fixed_ * X_PAR[n]';
    for(i in 1:PAR_size_)
      PAR_ranef[i] = 0;
    PAR_ranef = PAR_ranef + PAR_random_%[PAR_group_%[n]] * Z_PAR_%[n]';

    //PAR-link
    
    //likelihood
  }
}

model {
  for(i in 1:PAR_size)
    for(j in 1:X_PAR_ncol)
      PAR_fixed[i, j] ~ normal(PAR_prior_fixed_mu[i, j], PAR_prior_fixed_sd[i, j]);
  for(i in 1:PAR_size)
    for(j in 1:Z_PAR_ncol_%)
      PAR_sd_%[i, j] ~ cauchy(0, PAR_prior_scale_%[i, j]);
  L_corr_PAR_% ~ lkj_corr_cholesky(PAR_prior_nu_%);
  for(g in 1:PAR_group_max_%)
    PAR_z_%[g] ~ normal(0, 1);
  for(n in 1:N)
    counts[n] ~ multinomial(multinomial_p[n]);
}

generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));
}
