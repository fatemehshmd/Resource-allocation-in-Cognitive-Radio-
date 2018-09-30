%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get State and Psi
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang

% Update Date:
%   2011-02-18
%   2011-02-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, L, psi] = observe_state(BW, A, bitrate_set,beta_set,delta, coeff_psi) %, bitrate_set, beta_set, Psi_set, n_su
    % observed the updated states
    psi = zeros(size(A,1),1);
% 
    for ii=1:size(A,1)
        idx_beta = (beta_set==A(ii, 1));
        rate = (bitrate_set(idx_beta));
        psi(ii) = (1+BW./(rate.*(10.^(A(ii, 1)./10)))).^(-1);
%         psi(ii) = (1+1./(10.^(A(ii, 1)./10))).^(-1); %eq2
    end
    
%     psi = (1+1./(10.^(A./10))).^(-1); %eq2
%     psi = (1+BW./S(:, 1).*(10.^(S(:, 2)./10))).^(-1); %eq2
    psi_sum = sum(psi);
    if(psi_sum >= 1-delta)
        % the allocation failed, according to equation (3)
        I = 1;        
    else
        I = 0;        
    end
    
    %observe the interference to PU
    L = update_pu_sinr_level(coeff_psi, psi);
end