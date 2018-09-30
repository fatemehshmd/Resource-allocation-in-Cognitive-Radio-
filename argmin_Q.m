%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the q-table is stored in following way:
%
%  idx1 -- r
%  idx2 -- beta
%  idx3 -- I
%  idx4 -- L
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta, Q_max] = argmin_Q(Q_table, S,beta_set)
% S is the current state
% action is selected based on S

Q_sub(:,:) = Q_table(1,:, S(1)+1, S(2)+1);%by +1 we adjust the index

[Q_max, id_beta] = max(Q_sub); 
id_beta =find(Q_sub==Q_max);

% r=bitrate_set(id_r(id_beta));
beta = beta_set(id_beta(end));
end