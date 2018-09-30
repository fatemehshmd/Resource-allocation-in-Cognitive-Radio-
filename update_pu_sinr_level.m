%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hold position for state L, interference to PU
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = update_pu_sinr_level(coeff_psi, psi_su)
weighted_sum = sum(coeff_psi.* psi_su);
L = weighted_sum > 1;
end