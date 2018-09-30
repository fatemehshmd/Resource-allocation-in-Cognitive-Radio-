function [idx_beta, idx_I, idx_L] = find_idx(S, A,beta_set, I, L)

idx_beta = find(beta_set==A);

idx_I = find(I==S(1));
idx_L = find(L==S(2));

end