%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q-Learning algorithm for cooperative cognitive radio network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
flag_phase0 = 0;
flag_phase1 = 1; % add-one newcomer not learning 
flag_phase2= 1; % add-one newcomer learning from/using individual learning
flag_phase3 = 1; % add-one newcomer learning from/using individual learning dissimilar traffic
flag_phase4= 1; % add-one newcomer learning from/using individual learning similar traffic
flag_phase5= 1; % add-one newcomer learning from/using individual learning nearest neighbour
flag_phase6=1; %add-one newcomer learning from/using individual learning random neighbour
% system setupsDistortion_table
% Q-learning setups
learner_setup(1) = 0.01; %delta
learner_setup(2) = 0.8; %epsilon-greedy
learner_setup(3) = 0.4; %discounting factor
learner_setup(4) = 0.1; %learning rate
learner_setup(5) = 1; %cooperation flag
learner_setup (6)=0;


sys_setup(1) = 1e-9; %noise power 1nW
sys_setup(2) = 1e-2; %transmit power 10mW
sys_setup(3) = 10e6; %Hz
sys_setup(4) = 10^(1/10);  %primary user SINR requirement, from 1 dB

%create a set of possible transmit rate
[bitrate_set, beta_set, D_at_beta_br, Psi] = create_state_set(sys_setup(3));

Distortion_table{1} = bitrate_set;
Distortion_table{2} = beta_set;
Distortion_table{3} = D_at_beta_br;
Distortion_table{4} = Psi;

sz_I = 2; %binary set {0,1}
sz_L = 2; %binary set {0,1}
sz_br_set = size(bitrate_set, 2);
sz_beta_set = size(beta_set, 2);

flag_plot = 0;
t_max = 100;%maximum q-learning iterations
avg_time_max = 2000; %maximum q-learning iterations for a single setup
 
n_su =[4 8];  % No. of SUs

 learner_setup_new = learner_setup;
sz_n_su = size(n_su,2);

for i=1:sz_n_su
    su_traffic{i} = rand(1,n_su(i))>0.5;
end

avg_dis_exp = zeros(sz_n_su,1);
avg_conv_num_exp = zeros(sz_n_su,1);
avg_br_exp = zeros(sz_n_su,1);

avg_dis_new = zeros(sz_n_su,1);
avg_conv_num_new = zeros(sz_n_su,1);
avg_br_new = zeros(sz_n_su,1);

avg_dis_new_fi = zeros(sz_n_su,1);
avg_dis_new_fi_tr = zeros(sz_n_su,1);
avg_dis_new_fi_distr = zeros(sz_n_su,1);
avg_dis_new_fi_near = zeros(sz_n_su,1);
avg_conv_num_new_fi = zeros(sz_n_su,1);
avg_conv_num_new_fi_tr = zeros(sz_n_su,1);
avg_conv_num_new_fi_distr = zeros(sz_n_su,1);
avg_conv_num_new_fi_near = zeros(sz_n_su,1);
avg_br_new_fi_rnd = zeros(sz_n_su,1);
avg_br_new_fi = zeros(sz_n_su,1);
avg_br_new_fi_tr = zeros(sz_n_su,1);
avg_br_new_fi_distr = zeros(sz_n_su,1);
avg_br_new_fi_near = zeros(sz_n_su,1);
avg_br = zeros(sz_n_su,1);

avg_dis = zeros(sz_n_su,1);
avg_conv_num = zeros(sz_n_su,1);


avg_dis_PHYonly = zeros(sz_n_su,1);
avg_conv_num_PHYonly = zeros(sz_n_su,1);
avg_br_PHYonly = zeros(sz_n_su,1);

avg_dis_coop = zeros(sz_n_su,1);
avg_conv_num_coop = zeros(sz_n_su,1);
avg_br_coop = zeros(sz_n_su,1);

cong_rate = zeros(sz_n_su,1);
cong_rate_PHYonly = zeros(sz_n_su,1);
cong_rate_newcomer = zeros(sz_n_su,1);
cong_rate_newcomer_fi = zeros(sz_n_su,1);
cong_rate_newcomer_fi_tr = zeros(sz_n_su,1);
cong_rate_newcomer_fi_distr = zeros(sz_n_su,1);
cong_rate_newcomer_fi_near= zeros(sz_n_su,1);
cong_rate_coop = zeros(sz_n_su,1);
avg_dis_new_fi_rnd = zeros(sz_n_su,1);
avg_conv_num_new_fi_rnd = zeros(sz_n_su,1);
cong_rate_newcomer_fi_rnd = zeros(sz_n_su,1);
Gain_table = cell(3,1); %cell to store the gains randomly generated
Gain_table_new = cell(3,1); %cell to store the gains randomly generated with 1 addintional SU

%% Phase0
if flag_phase0
fprintf('individual learning - PHY Layer only\n');
learner_setup(5) = 0; %cooperation flag
for kk=1:sz_n_su    
    count = 0;
    sum_num_conv = 0;
    sum_dist = 0;
    sum_br = 0;
    for ii=1:avg_time_max        
        %generate positions and gains for the users
        [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));       
        Gain_table{1} = Gpu2pbs;
        Gain_table{2} = Gsu2pbs;
        Gain_table{3} = Gsu2sbs;        

        Q_table = zeros(n_su(kk), sz_br_set, sz_beta_set, sz_I, sz_L);
        expt = zeros(1, n_su(kk)); %accumulated expertness, has nothing to do with this algorithm      
        
        [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
            Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot);
        if(flag_converg == 1)
            %record the convergence
            sum_dist = sum_dist + avg_distortion;
            sum_num_conv = sum_num_conv + num_converge;
            sum_br = sum_br + avg_bitrate;
            count = count + 1;
        end
    end
    
    fprintf('Indivudual Learnning: No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop)\n', kk, n_su(kk), count);  
    
    if(count > 0)
        %avoid dividing by 0
        avg_dis_PHYonly(kk) = sum_dist / count;
        avg_conv_num_PHYonly(kk) = sum_num_conv / count;
        avg_br_PHYonly(kk) = sum_br / count;
        cong_rate_PHYonly(kk) = 1 - count / avg_time_max;
    end
end
end


%% Phase1
if flag_phase1
    %add-one new comer learning
    fprintf('new comer learning\n');
    vv = [];
    learner_setup(5) = 1;
    learner_setup_new(5) = 1;
    learner_setup_new (6) = 0;
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                %             expt_new(1:n_su(kk))=expt;
                expt_new = zeros(1, n_su(kk)+1); % fati
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1,sz_beta_set, sz_I, sz_L);
                
                
                learner_setup_new(6)=0;
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    sum_dist_new = sum_dist_new + avg_distortion;
                    vv = [vv num_converge];
                    sum_num_conv_new = sum_num_conv_new + (num_converge);
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop) , %d(doci)\n', kk, n_su(kk), count, count_new);
        
        if(count > 0)
            avg_dis_exp(kk) = sum_dist / count;
            avg_conv_num_exp(kk) = sum_num_conv / count;
            avg_br_exp(kk) = sum_bitrate / count;
        end
        
        if(count_new > 0)
            count_new;
            avg_dis_new(kk) = sum_dist_new / count_new;
            avg_conv_num_new(kk) = sum_num_conv_new / count_new;
            avg_br_new(kk) = sum_bitrate_new / count_new;
            cong_rate_newcomer(kk) = 1 - count_new / avg_time_max;
        end
       
    end
end 
%% Phase 2
if flag_phase2
    %add-one newcomer learning from/using individual learning
    fprintf('newcomer learning from/using individual learning\n');
    learner_setup(5) = 1;
    learner_setup_new(5) = 1; %cooperation flag
    learner_setup_new (6) = 1;
    vvv = [];
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new,p_su_x_new,p_su_y_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            
            
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                expt_new(1:n_su(kk))=expt;
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1, sz_beta_set, sz_I, sz_L);
                
                %keep the q-table of the old ones
                Q_table_new(1:n_su(kk),:, :, :) = Q_table;
                
                %copy the q-table into the new one's, notice that we are using
                %cooperative learning, all the q-table's are the same
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
                
                Q_table_new(n_su(kk)+1, :, :, :) = sum(Q_table,1)/n_su(kk);
                
               
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    vvv = [vvv num_converge];
                    sum_dist_new = sum_dist_new + avg_distortion;
                    sum_num_conv_new = sum_num_conv_new + num_converge;
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop), %d(doci)\n', kk, n_su(kk), count, count_new);
        
        
        if(count_new > 0)
            count_new;
            avg_dis_new_fi(kk) = sum_dist_new / count_new;
            avg_conv_num_new_fi(kk) = sum_num_conv_new / count_new;
            avg_br_new_fi(kk) = sum_bitrate_new / count;
            cong_rate_newcomer_fi(kk) = 1 - count_new / avg_time_max;
        end
    end
end
%% Phase 3
if flag_phase3
    % add-one newcomer learning from/using individual learning dissimilar traffic
    fprintf('newcomer learning from/using individual learning dissimilar traffic\n');
    learner_setup(5) = 1;
    learner_setup_new(5) = 1; %cooperation flag
    learner_setup_new (6) = 1;
    vvv = [];
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            
            
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                expt_new(1:n_su(kk))=expt;
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1, sz_beta_set, sz_I, sz_L);
                
                %keep the q-table of the old ones
                Q_table_new(1:n_su(kk),:, :, :) = Q_table;
                
                %copy the q-table into the new one's, notice that we are using
                %cooperative learning, all the q-table's are the same
                learner_setup_new (6) = 1;
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
             
                dissimilar_type_tr_users = find (su_traffic_type(1:end-1)~=su_traffic_type(end))';
                Q_table_new(n_su(kk)+1, :, :, :) = sum(Q_table(dissimilar_type_tr_users,:,:,:),1)/length(dissimilar_type_tr_users);
                
               
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    vvv = [vvv num_converge];
                    sum_dist_new = sum_dist_new + avg_distortion;
                    sum_num_conv_new = sum_num_conv_new + num_converge;
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop), %d(doci)\n', kk, n_su(kk), count, count_new);
        
        
        if(count_new > 0)
            count_new;
            avg_dis_new_fi_distr(kk) = sum_dist_new / count_new;
            avg_conv_num_new_fi_distr(kk) = sum_num_conv_new / count_new;
            avg_br_new_fi_distr(kk)= sum_bitrate_new / count;
            cong_rate_newcomer_fi_distr(kk) = 1 - count_new / avg_time_max;
        end
    end
end
%% Phase 4
if flag_phase4
    % add-one newcomer learning from/using individual learning similsr traffic
    fprintf('newcomer learning from/using individual learning similar traffic\n');
    learner_setup(5) = 1;
    learner_setup_new(5) = 1; %cooperation flag
    learner_setup_new (6) = 1;
    vvv = [];
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            
            
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                expt_new(1:n_su(kk))=expt;
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1, sz_beta_set, sz_I, sz_L);
                
                %keep the q-table of the old ones
                Q_table_new(1:n_su(kk),:, :, :) = Q_table;
                
                %copy the q-table into the new one's, notice that we are using
                %cooperative learning, all the q-table's are the same
                learner_setup_new (6) = 1;
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
             
                same_type_trf_users = find (su_traffic_type(1:end-1)==su_traffic_type(end))';
                Q_table_new(n_su(kk)+1, :, :, :) = sum(Q_table(same_type_trf_users,:,:,:),1)/length(same_type_trf_users);
                
               
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    vvv = [vvv num_converge];
                    sum_dist_new = sum_dist_new + avg_distortion;
                    sum_num_conv_new = sum_num_conv_new + num_converge;
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop), %d(doci)\n', kk, n_su(kk), count, count_new);
        
        
        if(count_new > 0)
            count_new;
            avg_dis_new_fi_tr(kk) = sum_dist_new / count_new;
            avg_conv_num_new_fi_tr(kk) = sum_num_conv_new / count_new;
            avg_br_new_fi_tr(kk)= sum_bitrate_new / count;
            cong_rate_newcomer_fi_tr(kk) = 1 - count_new / avg_time_max;
        end
    end
end
%% Phase 5
if flag_phase5
    % add-one newcomer learning from/using individual learning nearest
    % neighbour
    fprintf('newcomer learning from/using individual learning from nearest neighbour\n');
    learner_setup(5) = 1;
    learner_setup_new(5) = 1; %cooperation flag
    learner_setup_new (6) = 1;
    vvv = [];
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new,p_su_x_new,p_su_y_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            
            
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                expt_new(1:n_su(kk))=expt;
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1, sz_beta_set, sz_I, sz_L);
                
                %keep the q-table of the old ones
                Q_table_new(1:n_su(kk),:, :, :) = Q_table;
                
                %copy the q-table into the new one's, notice that we are using
                %cooperative learning, all the q-table's are the same
                learner_setup_new (6) = 1;
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
                
                su_num_nearest = findnearest(p_su_x, p_su_y,p_su_x_new,p_su_y_new);
             
                Q_table_new(n_su(kk)+1, :, :, :) = Q_table(su_num_nearest,:,:,:);
                
               
%                 same_type_trf_users = find (su_traffic_type(1:end-1)==su_traffic_type(end))';
%                 Q_table_new(n_su(kk)+1, :, :, :) = sum(Q_table(same_type_trf_users,:,:,:),1)/n_su(kk);
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    vvv = [vvv num_converge];
                    sum_dist_new = sum_dist_new + avg_distortion;
                    sum_num_conv_new = sum_num_conv_new + num_converge;
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop), %d(doci)\n', kk, n_su(kk), count, count_new);
        
        
        if(count_new > 0)
            count_new;
            avg_dis_new_fi_near(kk) = sum_dist_new / count_new;
            avg_conv_num_new_fi_near(kk) = sum_num_conv_new / count_new;
            avg_br_new_fi_near(kk)= sum_bitrate_new / count;
            cong_rate_newcomer_fi_near(kk) = 1 - count_new / avg_time_max;
        end
    end
end
%% Phase 6
if flag_phase6
    % add-one newcomer learning from/using individual learning random user
    % neighbour
    fprintf('newcomer learning from/using individual learning from random neighbour\n');
    learner_setup(5) = 1;
    learner_setup_new(5) = 1; %cooperation flag
    learner_setup_new (6) = 1;
    vvv = [];
    for kk=1:sz_n_su
        su_traffic_type = su_traffic{kk};%     su_traffic_type = [ones(1,n_su(kk)/2) zeros(1,n_su(kk)/2)];
        count = 0;
        sum_num_conv = 0;
        sum_dist = 0;
        sum_bitrate = 0;
        
        count_new = 0;
        sum_num_conv_new = 0;
        sum_dist_new = 0;
        sum_bitrate_new = 0;
        
        for ii=1:avg_time_max
            %generate positions and gains for the users
            [Gpu2pbs Gsu2pbs Gsu2sbs, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y] = randomize_G(n_su(kk));
            Gain_table{1} = Gpu2pbs;
            Gain_table{2} = Gsu2pbs;
            Gain_table{3} = Gsu2sbs;
            
            %generate positions by adding 1 SU
            [Gpu2pbs_new Gsu2pbs_new Gsu2sbs_new,p_su_x_new,p_su_y_new] = randomize_G2(n_su(kk)+1, p_su_x, p_su_y, p_pu_x, p_pu_y, p_sbs2pbs_x, p_sbs2pbs_y);
            Gain_table_new{1} = Gpu2pbs_new;
            Gain_table_new{2} = Gsu2pbs_new;
            Gain_table_new{3} = Gsu2sbs_new;
            
            Q_table = zeros(n_su(kk), sz_beta_set, sz_I, sz_L);
            expt = zeros(1, n_su(kk)); %accumulated expertness
            expt_new = zeros(1, n_su(kk)+1);
            
            
            %q-learning with cooperations, as the basis of the docitions
            [avg_distortion, num_converge, flag_converg, Q_table, expt, avg_bitrate] = q_learning_allocation2(n_su(kk), learner_setup, sys_setup, ...
                Q_table, expt, Gain_table, Distortion_table, t_max, flag_plot,su_traffic_type);
            
            if(flag_converg == 1)
                % only if it converges (no congestion occurs)
                sum_dist = sum_dist + avg_distortion;
                sum_num_conv = sum_num_conv + num_converge;
                sum_bitrate = sum_bitrate + avg_bitrate;
                count = count + 1;
                
                expt_new(1:n_su(kk))=expt;
                
                %add one SU to the system
                Q_table_new = zeros(n_su(kk)+1, sz_beta_set, sz_I, sz_L);
                
                %keep the q-table of the old ones
                Q_table_new(1:n_su(kk),:, :, :) = Q_table;
                
                %copy the q-table into the new one's, notice that we are using
                %cooperative learning, all the q-table's are the same
                learner_setup_new (6) = 1;
                su_traffic_type(n_su(kk)+1) = rand(1,1) > 0.5;
                
%                 su_num_nearest = findnearest_fati(p_su_x, p_su_y,p_su_x_new,p_su_y_new);
             su_num_nearest = randperm(n_su(kk),1);
                Q_table_new(n_su(kk)+1, :, :, :) = Q_table(su_num_nearest,:,:,:);
                
               
%                 same_type_trf_users = find (su_traffic_type(1:end-1)==su_traffic_type(end))';
%                 Q_table_new(n_su(kk)+1, :, :, :) = sum(Q_table(same_type_trf_users,:,:,:),1)/n_su(kk);
                [avg_distortion, num_converge, flag_converg_new, Q_table, expt_new, avg_bitrate] = q_learning_allocation2(n_su(kk)+1, learner_setup_new, ...
                    sys_setup, Q_table_new, expt_new, Gain_table_new, Distortion_table, t_max, flag_plot,su_traffic_type);
                
                if(flag_converg_new == 1)
                    vvv = [vvv num_converge];
                    sum_dist_new = sum_dist_new + avg_distortion;
                    sum_num_conv_new = sum_num_conv_new + num_converge;
                    sum_bitrate_new = sum_bitrate_new + avg_bitrate;
                    count_new = count_new + 1;
                end
            end
        end
        
        fprintf('No. of simulation: %d, No. of SUs: %d, Count of Convergence: %d(coop), %d(doci)\n', kk, n_su(kk), count, count_new);
        
        
        if(count_new > 0)
            count_new;
            avg_dis_new_fi_rnd(kk) = sum_dist_new / count_new;
            avg_conv_num_new_fi_rnd(kk) = sum_num_conv_new / count_new;
            avg_br_new_fi_rnd(kk)= sum_bitrate_new / count;
            cong_rate_newcomer_fi_rnd(kk) = 1 - count_new / avg_time_max;
        end
    end
end
%% PLOT
n_su_new = n_su+1;

figure;
plot(n_su_new, avg_dis_new, 'o-','linewidth',3,'MarkerSize',12); hold on;
% plot(n_su, avg_dis_coop, '*-'); hold on;
plot(n_su_new, avg_dis_new_fi, 'p-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_dis_new_fi_tr, 'co-','linewidth',3,'MarkerSize',12); hold on;
plot(n_su_new, avg_dis_new_fi_distr, 'r*-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_dis_new_fi_near, 'r+-','linewidth',3,'MarkerSize',12);
plot(n_su_new, avg_dis_new_fi_rnd, '*-','linewidth',3,'MarkerSize',12); 
grid;
ylabel('Average MOS');
xlabel('Number of Secondary Users');
legend( 'Newcomer f/individual','Newcomer f/individual avr','Newcomer f/individual similar traffic','Newcomer f/individual dissimilar traffic','Newcomer f/individual nearest neighbour','Newcomer f/individual random neighbour');
set(gca,'FontSize',20,'fontweight','bold');

figure;
plot(n_su_new, avg_br_new, 'o-','linewidth',3,'MarkerSize',12); hold on;
% plot(n_su, avg_br_coop, '*-'); hold on;
plot(n_su_new, avg_br_new_fi, 'p-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_br_new_fi_tr, 'c+-','linewidth',3,'MarkerSize',12);
plot(n_su_new, avg_br_new_fi_distr, 'r*-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_br_new_fi_near, 'ro-','linewidth',3,'MarkerSize',12); hold on;
plot(n_su_new, avg_br_new_fi_rnd, '*-','linewidth',3,'MarkerSize',12);
grid;
ylabel('Average Bitrate');
xlabel('Number of Secondary Users');
% legend('Individual Learning','Cooperative Learning', 'Cooperative Initialization for New Comer', 'New Comer Learning', 'Newcomer f/individual','PHY only individual');
legend('Newcomer f/individual', 'Newcomer f/individual avr','Newcomer f/individual similar traffic','Newcomer f/individual dissimilar traffic','Newcomer f/individual nearest neighbour','Newcomer f/individual random neighbour');
set(gca,'FontSize',20,'fontweight','bold');

figure;
plot(n_su_new, avg_conv_num_new, 'o-','linewidth',3,'MarkerSize',12); hold on;
% plot(n_su, avg_conv_num_exp, 'ro-'); hold on
plot(n_su_new, avg_conv_num_new_fi, '*-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_conv_num_new_fi_tr, 'c+-','linewidth',3,'MarkerSize',12);
plot(n_su_new, avg_conv_num_new_fi_distr, 'rp-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_conv_num_new_fi_near, '*-','linewidth',3,'MarkerSize',12);hold on;
plot(n_su_new, avg_conv_num_new_fi_rnd, '*-','linewidth',3,'MarkerSize',12)
grid;
ylabel('Average Iteration Number of Convergence');
xlabel('Number of Secondary Users');
% legend('Individual Learning','Cooperative Learning', 'Cooperative Initialization for New Comer', 'New Comer Learning', 'Newcomer f/individual','PHY only individual');
legend('Newcomer f/individual','Newcomer f/individual avr','Newcomer f/individual similar traffic','Newcomer f/individual dissimilar traffic','Newcomer f/individual nearest neighbour','Newcomer f/individual random neighbour');
set(gca,'FontSize',20,'fontweight','bold');

figure;
% plot(n_su, cong_rate, 'o-',n_su, cong_rate_coop, '*-',n_su_new, cong_rate_newcomer, 's-',n_su_new, cong_rate_newcomer_fi, 'p-'); grid;
plot(n_su_new,cong_rate_newcomer, 'o-',n_su_new, cong_rate_newcomer_fi, 'p-',n_su_new,cong_rate_newcomer_fi_tr, '*-',n_su_new, cong_rate_newcomer_fi_distr, 's-',n_su_new, cong_rate_newcomer_fi_near, 'o-',n_su_new, cong_rate_newcomer_fi_rnd,'c-*','linewidth',3,'MarkerSize',12); grid;

ylabel('Congestion Rate');
xlabel('Number of Secondary Users');
% legend('Individual Learning', 'Cooperative Learning', 'Newcomer Learning', 'Newcomer f/individual','PHY only individual');
legend('Newcomer f/individual','Newcomer f/individual avr','Newcomer f/individual similar traffic','Newcomer f/individual dissimilar traffic','Newcomer f/individual nearest neighbour','Newcomer f/individual random neighbour');
set(gca,'FontSize',20,'fontweight','bold');

