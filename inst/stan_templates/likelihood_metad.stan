// meta-d' likelihood
shift = 0.5 * stim_sign[n]; // when stim == 1 shift == -0.5
normalization[1] = Phi(gamma_[Kb2] - shift * delta_[1]) / Phi(gamma_[Kb2] - shift * delta_[2]);
normalization[2] = Phi(-(gamma_[Kb2] - shift * delta_[1])) / Phi(-(gamma_[Kb2] - shift * delta_[2]));
multinomial_p[n, 1] = Phi(gamma_[1] - shift * delta_[2]) * normalization[1];
for(k in 2:(K - 1))
  if(k < (Kb2 + 1)){
    multinomial_p[n, k] = (Phi(gamma_[k] - shift * delta_[2]) - Phi(gamma_[k - 1] - shift * delta_[2])) * normalization[1];
  }else{
    multinomial_p[n, k] = (Phi(gamma_[k] - shift * delta_[2]) - Phi(gamma_[k - 1] - shift * delta_[2])) * normalization[2];
  }
multinomial_p[n, K] = Phi(-(gamma_[K - 1] - shift * delta_[2])) * normalization[2];
// meta-d' likelihood end
