//gamma_link:log_ratio
gamma_[Kb2] = gamma[Kb2]; // main criterion, i.e., bias
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
