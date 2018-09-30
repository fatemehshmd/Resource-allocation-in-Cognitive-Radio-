%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q-Learning algorithm for cooperative source-channel flow control
% Author:
%   Wenbo Wang, Andres, Kwasinski
% Time of Creation:
%   2010-02-13
% Time of Modification:
%   2010-02-15
%   2010-02-18
%   2010-02-20
%   2010-02-23
%   2010-02-24
%   2010-02-25
%   2010-02-27
%`  2010-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg_distortion, num_converge, flag_converg, Q_table_ret, expt_ret, avg_bitrate] = q_learning_allocation2(n_su, learner_setup, sys_setup, ...
    Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type)

delta = learner_setup(1); %ceiling for SU SINR

% Q-learning setups
epsilon = learner_setup(2); %epsilon-greedy
gamma = learner_setup(3); %discounting factor
alpha = learner_setup(4); %learning rate
flag_coop = learner_setup(5); %cooperation flag

lambda_rate = 1;
lambda_decay = 0.9;

%create a set of possible transmit rate
%[bitrate_set, beta_set, D_at_beta_br] = create_state_set();

sigma2 = sys_setup(1);
P0 = sys_setup(2); %transmit power
BW = sys_setup(3); %Hz
beta_0 = sys_setup(4); %from 13 dB

Gpu2pbs = Gain_table{1};
Gsu2pbs = Gain_table{2};
Gsu2sbs = Gain_table{3};
%Initialize the channel gain with respect to positions
%[Gpu2pbs Gsu2pbs Gsu2sbs] = randomize_G(n_su);

coeff_psi = Gsu2pbs.*(sigma2+Gpu2pbs*P0)./(Gsu2sbs.*(Gpu2pbs.*P0./beta_0-sigma2))+1;

bitrate_set = Distortion_table{1};
beta_set = Distortion_table{2};
D_at_beta_br = Distortion_table{3};
Psi_set = Distortion_table{4};



%Initialization before iterations
% expt = zeros(1, n_su); %accumulated expertness
c = zeros(1, n_su); % immediate cost function
weight_exp = zeros(1, n_su);
D = zeros(1, n_su); % part of the cost function (distortion)

% Initialize the state-actions, {r, beta, I, L}, I={0,1}, L={0,1};
S = zeros(n_su, 2);
S_new = zeros(n_su, 2);
% A = repmat([beta_set(end)], n_su, 1); %old actions
A_new = zeros(n_su, 1); %new actions
Current_Q = zeros(1, n_su);


 % Set the required transmit rate randomly
[tran_rate_level, beta_idx] = set_reqired_trans_rate(bitrate_set, n_su);
A = ([beta_set(beta_idx)]'); %old actions
A = repmat([beta_set(end)], n_su, 1); %old actions
    
    

if  learner_setup(6)
    for ii=1:n_su
        
        [A(ii),Q_max] = argmin_Q(Q_table(ii,:,:,:), S(ii,:),beta_set);%find an action
    end
end


%fill up Q-table with the first observation
[I, L, psi_su] = observe_state(BW, A ,bitrate_set,beta_set, delta, coeff_psi); %bitrate_set, beta_set, Psi_set, n_su (obsolete parameters)
S(:,1) = I;
S(:,2) = L;

% for testing purpose TO BE DELETED
% beta_idx(:,:)=1;

% Q-learning begin

% Building up Q-table hierarchies for each SU, permutations of state and
% action (initial value = 0), dementionality is:
% #(bit_rate)*#(beta_vector)*#(I)*#(L)
%Q_table = zeros(n_su, sz_br_set, sz_beta_set, sz_I, sz_L);

%set Q_tables elements lower than the minimum transmit rate level to some
%large value
% if ~ learner_setup(6)
%     if I+L>=1
%         for ii=1:n_su
%             if(beta_idx(ii)>1)
%                 Q_table(ii, 1:beta_idx(ii)-1, :, :) = -1;
%             end
%         end
%     end
% end
if ~ learner_setup(6)
    
    for ii=1:n_su
        D(ii) = get_distortion(A(ii), beta_set, D_at_beta_br,su_traffic_type(ii));
        %     c(ii) = D(ii) + lambda(ii)*psi_su(ii) + (I+L)*(D_at_beta_br(1,1)+lambda(ii)*Psi_set(id_end_beta));
        c(ii) = D(ii);
        
        [idx_beta, idx_I, idx_L] = find_idx(S(ii,:), A(ii), beta_set, [0,1], [0,1]);
        
        if(Q_table(ii, idx_beta, idx_I, idx_L) == 0)
            if I+L<1
                Q_table(ii, idx_beta, idx_I, idx_L) = c(ii);
            else
                Q_table(ii, idx_beta, idx_I, idx_L) = - 0.03;
            end
        end
        
    end
end


lambda1 = 1;%50+(100-50).*rand(1, 1); %randomize lambda1
lambda2 = 1;%50+(100-50).*rand(1, 1); %randomize lambda2
lambda = lambda1 + coeff_psi.*lambda2; %

%sum of distortion
sum_D = zeros(1, t_max);
sum_IL = zeros(1, t_max);



id_end_r = size(bitrate_set, 2);
id_end_beta = size(beta_set, 2);


flag_same_state = 0;
count_same_state = 0;

state_trend =[];
actions = [];

%% update Q-values
for tt=1:t_max
    
    Current_Q(:,:) = 0;
    D(:,:) = 0;
    c(:,:) = 0;
    
    %update espilon
    epsilon = 1-(1-epsilon)/tt; % Question
    eps = rand(n_su, 1); %to determine if each of the SUs enter the learning process
    indicator = eps<=epsilon;
    
    for ii=1:n_su
        if(indicator(ii) == 1)
            %learn from Q value, to find the new action (r_i, beta_i) which
            %minimizes the old Q-value table
            tt;
            [A_new(ii),Q_max] = argmin_Q(Q_table(ii,:,:,:), S(ii,:),beta_set);%find an action
        else
            %            idx_br_i = randi([1,size(bitrate_set,2)], 1);
            idx_beta_i = randi([1,size(beta_set,2)], 1);
            %            A_new(ii,1) = bitrate_set(idx_br_i);
            A_new(ii) = beta_set(idx_beta_i);
            %A_new(ii,:) = A(ii,:);
        end
        
    end
    
    
    % observed the updated states
    [I, L, psi_su] = observe_state(BW, A_new,bitrate_set,beta_set, delta, coeff_psi);
    
    S_new(:,1) = I;
    S_new(:,2) = L;
    
    
    
    % get immediate cost value c by observation
    for ii=1:n_su
        %if(indicator(ii) == 1)
        D(ii) = get_distortion(A_new(ii,1), beta_set, D_at_beta_br,su_traffic_type(ii));
        
        %get immediate c from equation (12) with r and beta
        if(I+L==0)
            c(ii) = D(ii);
        else
            c(ii) = -5; % Question
            %             c(ii) = - 0.5;
        end
        %  c(ii) = D(ii) + lambda(ii)*psi_su(ii) + (I+L)*(D_at_beta_br(1,1)+lambda(ii)*Psi_set(id_end_r,id_end_beta));  %
        
        %update q-table
        %find the minimum Q-value for new state
        [ beta, Q_min] = argmin_Q(Q_table(ii,:,:,:), S_new(ii,:), beta_set);%retrieve the q-values for updating
        
        %combination of old state and new actions
        [idx_beta, idx_I, idx_L] = find_idx(S(ii,:), A_new(ii,:), beta_set, [0,1], [0,1]);
        
        %update with function (9)
        Current_Q(ii) = (1-alpha)*Q_table(ii, idx_beta, idx_I, idx_L) + alpha*(c(ii)+gamma*Q_min);
        Q_table(ii,idx_beta, idx_I, idx_L) = Current_Q(ii);
        %end
    end
    
    
    
    
    %update expertness
    expt = expt + 1./(lambda'.*c);%gamma*expt+c;  % what is that
    
    if(flag_coop==1)
        expt_sum = sum(expt);
        for ii=1:n_su
            weight_exp(ii) = expt(ii)/expt_sum;
        end
        
        % cooperation with other SU
        for ii=1:n_su
            [idx_beta, idx_I, idx_L] = find_idx(S(ii,:), A_new(ii,:), beta_set, [0,1], [0,1]);
            
            %that should be weighed q-table cell of all the other's
            Q_table(ii,idx_beta, idx_I, idx_L) = weight_exp*Q_table(:, idx_beta, idx_I, idx_L);
        end
    end
    
    
    
    
    %update lamdbda towards the direction of gradient
    %     if (I+L == 0 ) %&& tt >= size(Distortion_table{2},1)*size(Distortion_table{2},2)
    %         % One way of updating the results
    %         lambda1 = lambda1 + lambda_rate * sum(psi_su.^2)/n_su;
    %         lambda2 = lambda2 + lambda_rate * sum((psi_su.*coeff_psi).^2) / n_su;
    %
    %         % Another way of updating the results
    %         %lambda1 = lambda1 + lambda_rate * sqrt(sum(psi_su.^2));
    %         %lambda2 = lambda2 + lambda_rate * sqrt(sum((psi_su.*coeff_psi).^2));
    %
    %         lambda = lambda1 + coeff_psi.*lambda2;
    %         lambda_rate = lambda_rate*lambda_decay;
    %     end
    
    
    
    
    
    if(sum(sum(S==S_new))==2*n_su && sum(sum(A == A_new)) == n_su)
        flag_same_state = 1;
    else
        flag_same_state = 0;
    end
    
    
    
    if(flag_same_state == 1)
        count_same_state = count_same_state + 1;
    else
        count_same_state = 0;
    end
    
    %update state
    S = S_new;
    A = A_new;
    actions(tt,1:n_su) = [A_new']; % first col user one, sec col user two
    state_trend (tt,1:n_su) = sum(S_new');
    %record sum of distortion
    sum_D(tt) = sum(D);
    sum_IL(tt) = I+L; %  Question why
    D = [];
    
    if(I+L == 0 && count_same_state >= 4) % Question Condition for convergence???
        break; % Question
    end
    
end

Q_table_ret = Q_table;
expt_ret = expt;

if(I+L == 0 && count_same_state >= 4 )
    flag_converg = 1;
    rate = zeros(size(n_su));
    for i=1:n_su
        rate(i)= bitrate_set(beta_set==A(i));
        MOS (i) = D_at_beta_br(beta_set==A(i));
    end
    avg_bitrate = sum(rate) / n_su;
%     avg_distortion = sum (MOS) /n_su;
    avg_distortion = sum_D(tt)/n_su;
    num_converge = tt;
    A;
else
    flag_converg = 0;
    avg_bitrate = 0;
    avg_distortion = 0;
    num_converge = 0;
end

if(flag_plot == 1)
    figure
    plot(1:tt, actions(:,2),'-*');
    xlabel('Iteration');
    ylabel('Action');
end

end
