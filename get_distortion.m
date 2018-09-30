%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The distortion is obtained with equation 14.30, chapter 14.6.3 in book 
% "Cooperative Communications and Networking" by K.J. Liu, A.K. Sadek,
% W. Su and A. Kwasinski
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = get_distortion(beta,beta_set, D_at_beta_br,p)
%     idx_r = find(bitrate_set==r);
%     if(idx_r == 0)
%        disp('incorrect index'); 
%     end
    
    idx_beta = find(beta_set==beta);
    if(idx_beta == 0)
       disp('incorrect index'); 
    end
    
    D = D_at_beta_br(1,idx_beta,p+1);
end