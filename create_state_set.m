%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The distortion is obtained with equation 14.30, chapter 14.6.3 in book 
% "Cooperative Communications and Networking" by K.J. Liu, A.K. Sadek,
% W. Su and A. Kwasinski
%
% Author: A. Kwasisnki
% Last Update by: Wenbo Wang
%
% Last Update: 2011-02-26
%              2011-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [br_set,Beta_dB, D_at_beta_br, Psi] = create_state_set(BW)

% RCPC setting 
NP = 100; 
a = 6.0769;
BER = 10^(-5);
% BER = 0.2*exp(-1.5.*Beta_dB) % appendix A of paper performance of adaptive mod in cellular systems
k = 1.5/(-log(5*BER));
PER = 1- (1-BER).^ NP;

Psi = [];

% channel_code_rates= [8/9,4/5,2/3,4/7,1/2,4/9,4/10,4/11,1/3,4/13,2/7,4/15,1/4];
% br_set = [0.09765625e5, 0.1953125e5, 0.390625e5, 0.78125e5, 1.5625e5, 3.125e5];

% channel_code_rates =[4/10,4/11,1/3,4/13,2/7,4/15,1/4];
% br_set = [0.78125e5, 1.5625e5, 3.125e5];

% Beta_dB = zeros(length(br_set), length(channel_code_rates));
% Beta = zeros(1, length(channel_code_rates)); 

% Psi = zeros(length(br_set), length(channel_code_rates));
%%
Beta_dB = [-5:2:15];
Beta = 10.^(Beta_dB./10);
br_set = BW * log2(1+ (k.* Beta))./1e3;
D_at_beta = zeros(1,length(br_set),2);      
% MOS data
D_at_beta_br(1,:,1) = 1.4 .* log10( 0.4780 .* br_set); % updated version of paper A survey on parametric QoE estimation for popular services
% MOS video
br_set = BW * log2(1+(k.* Beta));
PSNR = 10.4 * log10 (br_set)- 28.7221;
% Beta_dB = 10.4 * log10(br_set) - 28.7221;
% D_at_beta_br(1,:,2) = p(2)./(1 + exp(-p(3)*(PSNRperVidPiTr-p(1))));
D_at_beta_br(1,:,2) = 6.6431./(1+ exp(-0.1344 * (PSNR - 30.4264)));
%  D_at_beta_br(1,:,2) = D_at_beta_br(1,:,1);
           
        
% PER (1,:) = 1- (1- (a/2)* erfc(sqrt(30 * exp(-3 .* channel_code_rates).* Beta))).^ (NP * channel_code_rates * br_set(1));

% Psi(ii,jj) = (1+BW./(br_set(ii).*(10.^(Beta_dB(ii,jj)./10)))).^(-1);
% Psi = (1+1./((10.^(Beta_dB./10)))).^(-1);
Psi = (1+BW./(br_set.*(10.^(Beta_dB./10)))).^(-1);

%% second alternative
% br_set = [0.4:0.4:3]*1e6;
% Beta = (2.^(br_set./BW)-1)./k;
% Beta_dB = 10 * log10(Beta);
% D_at_beta = zeros(1,length(br_set),2);      
% % MOS data
% D_at_beta_br(1,:,1) = 1.3619 .* log10( 0.6780 .* br_set); % updated version of paper A survey on parametric QoE estimation for popular services
% % MOS video
% % Beta_dB = 10.4 * log10(br_set) - 28.7221;
% % D_at_beta_br(1,:,2) = 6.6431./(1+ exp(-0.1344 * (Beta_dB - 30.4264)));
% D_at_beta_br(1,:,2) = D_at_beta_br(1,:,1);
% 
% Psi = (1+BW./(br_set.*(10.^(Beta_dB./10)))).^(-1);
end 