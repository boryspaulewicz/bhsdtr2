// ordinal likelihood
for(k in 1:(K - 1))
  multinomial_cum[k + 1] = Phi(gamma_[k] - eta_[1]);
multinomial_cum[1] = 0;
multinomial_cum[K + 1] = 1;
for(k in 1:K)
  multinomial_p[n, k] = multinomial_cum[k + 1] - multinomial_cum[k];
// ordinal likelihood end
