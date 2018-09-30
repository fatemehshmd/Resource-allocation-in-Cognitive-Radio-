function [trans_rate_level, idx_in_br_set] = set_reqired_trans_rate(bitrate_set, n_su)

% if we set them randomly
sz_br_set = size(bitrate_set, 2);
idx_in_br_set = randi([3,sz_br_set], 1, n_su);
trans_rate_level = bitrate_set(idx_in_br_set);
end