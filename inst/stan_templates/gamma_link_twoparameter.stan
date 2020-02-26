//gamma_link:twoparameter
gamma_[Kb2] = gamma[1]; // main criterion i.e. bias
if(K > 2){
  for(k in 1:(Kb2 - 1)){
    gamma_[Kb2 + k] = gamma_[Kb2 + k - 1] + exp(gamma[2]);
    gamma_[Kb2 - k] = gamma_[Kb2 - k + 1] - exp(gamma[2]);
  }
 }
