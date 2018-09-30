%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geographical Distribution Generator for SU, SBS, PU, PBS
%
% Author: Wenbo Wang
% Last Update by: Wenbo Wang

% Update Date:
%   2011-02-18
%   2011-02-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gpu2pbs, Gsu2pbs, Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(Nsu)
n = 2.8;%large scale fading

d0 = 10;%reference distance

% Suppose PBS position is (0,0), randomize the SBS position, SU and PU
% position respectively

distance_near = -200;
distance_far = 200;
p_pu_x = distance_near+(distance_far-distance_near).*rand(1, 1);%coordinate of PU
p_pu_y = distance_near+(distance_far-distance_near).*rand(1, 1);%coordinate of PU

%SBS is within a circle of radius 200
radius = 2000;% fati
p_sbs2pbs_x = -radius+(radius+radius)*rand(1, 1);%coordinate of SU
p_sbs2pbs_y = -radius+(radius+radius)*rand(1, 1);%coordinate of SU
% p_sbs2pbs_x = -radius+(radius+radius)*0.8;%coordinate of SU
% p_sbs2pbs_y = -radius+(radius+radius)*0.8;%coordinate of SU
%SUs are distributed around the SBS
% distance_near_x = p_sbs2pbs_x-250;
% distance_far_x = p_sbs2pbs_x+250;
% distance_near_y = p_sbs2pbs_y-250;
% distance_far_y = p_sbs2pbs_y+250;
r=1000;
distance_near_x = p_sbs2pbs_x-r;
distance_far_x = p_sbs2pbs_x+r;
distance_near_y = p_sbs2pbs_y-r;
distance_far_y = p_sbs2pbs_y+r;

p_su_x = distance_near_x+(distance_far_x-distance_near_x).*rand(Nsu, 1);%coordinate of SU
p_su_y = distance_near_y+(distance_far_y-distance_near_y).*rand(Nsu, 1);%coordinate of SU

%Gains for each transeiver pair
Gpu2pbs = 1./((sqrt(p_pu_x.^2+p_pu_y.^2)./d0).^n); %channel gain -pu to pbs
Gsu2pbs = 1./((sqrt(p_su_x.^2+p_su_y.^2)./d0).^n); %channel gain -su to pbs

Gsu2sbs = 1./((sqrt((p_su_x-p_sbs2pbs_x).^2+(p_su_y-p_sbs2pbs_y).^2)./d0).^n); %channel gain -su to sbs


% figure;
% plot(p_su_x, p_su_y, 'ro', 'MarkerSize', 6); hold on;
% plot(p_pu_x, p_pu_y, 'bo', 'MarkerSize', 6); hold on;
% plot(0, 0, 'b*', 'MarkerSize', 6); hold on;
% plot(p_sbs2pbs_x, p_sbs2pbs_y, 'r*', 'MarkerSize', 6); hold on;
% title('Position Distribution of the SUs and PU');
% legend('SU','PU','PBS','SBS');

end