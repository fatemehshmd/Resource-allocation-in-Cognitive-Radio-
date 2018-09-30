%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geographical Distribution Generator for SU, SBS, PU, PBS
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang

% Update Date:
%   2011-02-18
%   2011-02-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gpu2pbs Gsu2pbs Gsu2sbs,p_su_x_new,p_su_y_new] = randomize_G2(Nsu, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y)
n = 2.8;%large scale fading
d0 = 10;%reference distance

% Suppose PBS position is (0,0), randomize the SBS position, SU and PU
% position respectively

%SUs are distributed around the SBS
r= 1000;
distance_near_x = p_sbs2pbs_x-r;
distance_far_x = p_sbs2pbs_x+r;
distance_near_y = p_sbs2pbs_y-r;
distance_far_y = p_sbs2pbs_y+r;

p_su_x_new = distance_near_x+(distance_far_x-distance_near_x).*rand(Nsu-size(p_su_x,1), 1);%coordinate of SU
p_su_y_new = distance_near_y+(distance_far_y-distance_near_y).*rand(Nsu-size(p_su_x,1), 1);%coordinate of SU

p_su_x = [p_su_x; p_su_x_new];
p_su_y = [p_su_y; p_su_y_new];

%Gains for each transeiver pair
Gpu2pbs = 1./((sqrt(p_pu_x.^2+p_pu_y.^2)./d0).^n); %channel gain -pu to pbs
Gsu2pbs = 1./((sqrt(p_su_x.^2+p_su_y.^2)./d0).^n); %channel gain -su to pbs

Gsu2sbs = 1./((sqrt((p_su_x-p_sbs2pbs_x).^2+(p_su_y-p_sbs2pbs_y).^2)./d0).^n); %channel gain -su to sbs

end