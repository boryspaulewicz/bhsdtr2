//gamma_link:softmax
gamma_ = thresholds_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma, 0))), K - 1));
