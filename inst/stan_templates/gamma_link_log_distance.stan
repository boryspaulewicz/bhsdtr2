//gamma_link:log_distance
gamma_[Kb2] = gamma[Kb2]; // main criterion i.e. bias
if(K > 2){
  for(k in 1:(Kb2 - 1))
    gamma_[Kb2 - k] = gamma_[Kb2 - k + 1] - exp(gamma[Kb2 - k]);
  for(k in (Kb2 + 1):(K - 1))
    gamma_[k] = gamma_[k - 1] + exp(gamma[k]);
 }

