clear

C_prev = 0;
C = 200; % total number of cycles
check_cycle = C;
test_rep_num = 3;
spike_frac = 0; % fraction of H pure culture spiked in
spike_test = [0.35 0.7];
sl = length(spike_test);
% upper bound for gH_max. for gH_max_Bound = 0.3, set V = 1.
gH_max_Bound = 0.8;
% the factor for amount of R(0)
V = 10;
% minimal number of Adults allowed to reproduce. comm_type_num = 1 for the
% top-dog strategy, comm_type_num = n for the top n% strategy.
comm_type_num = 2;
% mutation rate corresponding to effective mutation rate of 2e-3
% to turn off the mutation, set mut_rate=0.
mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
Pn_sig = 100; % std of the noise 40 correpondes to 5% and 80 corresponds to 10%
% reproducing method
% @pipette_SP: reproduce through pipetting, allowing both BM(0) and
%          phi_M(0) to fluctuate
% @cell_sort_SP: reproduce through cell-sorting, so that
%            the biomass in the Newborns are fixed BM(0)=BM_target and
%            phi_M(0)=phi_M(T) of the parent Adults
repro_method = @pipette_spike_SP;

N = 100; % number of communities within a cycle
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
T0 = 17;
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4*V; % maximal number of cells in the community
nb_popul = BM_target*2;
t_bin = 0.05; % time step in the simulation
pcs=1e-15; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

% ancestral or evolved parameters shown in Table 1
gM_max_start = 0.7/1.2; % max growth rate of M
gH_max_start = 0.3/1.2; % max growth rate of H
dM = 3.5e-3; % death rate of M
dH = 1.5e-3; % death rate of H
fp_start = 0.1; % fp at the beginning of selection
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR_start = 1 ;
K_HR_start = 1 ;
K_MB_start = 500/3 ;

% evolutionary bounds of the phenotypes
gM_max_Bound = 0.7;
fp_Bound = 1; % fp is between 0 and 1
K_MB_Bound = 100/3 ;
K_MR_Bound = 0.33 ;
K_HR_Bound = 0.2 ;


% R(0)=1 for each V
R0 = V;
% if K_MR, K_HR, K_MB is larger than K_singular value, they are effectively
% infinity
K_singular = 1e3;


rng('shuffle');
comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'gM_max', zeros(max_popul,1), 'K_MB', zeros(max_popul,1), 'K_MR', zeros(max_popul,1),...
    'gH_max', zeros(max_popul,1), 'K_HR', zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
newborn_struct = struct('M_L',zeros(nb_popul,1),'H_L',zeros(nb_popul,1),'fp',zeros(nb_popul,1),...
    'gM_max', zeros(nb_popul,1), 'K_MB', zeros(nb_popul,1), 'K_MR', zeros(nb_popul,1),...
    'gH_max', zeros(nb_popul,1), 'K_HR', zeros(nb_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target, ...
    'gH_max_start', gH_max_start, 'K_HR_start', K_HR_start, 't_bin', t_bin,...
    'gM_max_Bound', gM_max_Bound, 'gH_max_Bound', gH_max_Bound,...
    'fp_Bound', fp_Bound, 'K_MB_Bound', K_MB_Bound,'K_MR_Bound', K_MR_Bound, ...
    'K_HR_Bound', K_HR_Bound, 'R0', R0, 'K_singular', K_singular, 'dM', dM, 'dH', dH, ...
    'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH, 'mut_rate', mut_rate);
%%
if C_prev > 0
    load('comm_all/newborns')
else
    fp_init = zeros(nb_popul, 1);
    gM_max_init = zeros(nb_popul, 1);
    K_MB_init = zeros(nb_popul, 1);
    K_MR_init = zeros(nb_popul, 1);
    gH_max_init = zeros(nb_popul, 1);
    M_L_init = zeros(nb_popul, 1);
    K_HR_init = zeros(nb_popul, 1);
    H_L_init=zeros(nb_popul, 1);
    % death_probability=zeros(max_popul,1);
    M_counter = 60 ;
    H_counter = BM_target - M_counter;
    M_L_init(1 : M_counter) = 1;
    H_L_init(1 : H_counter) = 1;
    fp_init(1 : M_counter) = fp_start;
    gM_max_init(1 : M_counter) = gM_max_start;
    K_MB_init(1 : M_counter) = K_MB_start;
    K_MR_init(1 : M_counter) = K_MR_start;
    gH_max_init(1 : H_counter)=gH_max_start;
    K_HR_init(1 : H_counter) = K_HR_start;
    
    newborns(1, N)...
        =struct('M_L', M_L_init,'H_L', H_L_init,'fp', fp_init,...
        'gM_max', gM_max_init, 'K_MB', K_MB_init, 'K_MR', K_MR_init,...
        'gH_max', gH_max_init, 'K_HR', K_HR_init,...
        'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
        'P',0,'parentnum',0,'rseed',uint32(0));
    rseed=randi(2^32-1,N,1,'uint32');
    for ri=1:N
        newborns(ri).rseed=rseed(ri);
    end
end

n = C_prev+1;
% for n = C_prev+1 : max(C, check_cycle)
while n <= max(C, check_cycle)
    % create a folder Cn to save the results of the nth cycle
    folder_name1 = ['C' num2str(n)];
    if ~exist(folder_name1, 'dir')
        mkdir(folder_name1)
    end
    folder_name2=['C' num2str(n) '/comm_all'];
    if ~exist(folder_name2, 'dir')
        mkdir(folder_name2)
    end
    save([folder_name2 '/newborns'],'newborns');
    comm_all(1,1:N) = comm_struct;
    % rep is the index of communities within one cycle
    parfor rep = 1:N
        %   for rep = 1:2
        comm_all(rep) = simu_one_comm(newborns(rep), comm_struct, const_struct);
    end
    distrng=rng;
    % add stochastic noise to the community function
    Pn = normrnd(0, Pn_sig, size([comm_all.P]));
    P_all = [comm_all.P];
    [~, I] = sort([comm_all.P] + Pn, 'descend');
    % I=randperm(comm_type_num*comm_rep_num);
    check_flag = false;
    if abs(n - check_cycle + 2) < 0.1
        sel_counter = 0;
        rep_counter = 0;
        % rep_num_M is a matrix for the number of newborns for selection
        % and test for each chosen adult
        rep_num_M = zeros(N, sl+1);
        for i = 1 : comm_type_num*comm_rep_num
            if rep_counter >= comm_type_num*comm_rep_num
                break
            end
            BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
            % number of test comms the Adult can contribute
            test_rep_temp = floor((BM - comm_rep_num*BM_target*(1-spike_frac))...
                    /BM_target/sum(1-spike_test));
            % if the number of test comms is too small, increase the test cycle by 1           
            if test_rep_temp < round(comm_rep_num/5)
                check_cycle = check_cycle+1;
                check_flag = false;
                break
            else
                rep_num_M(i, 1) = comm_rep_num;
                rep_num_M(i, 2:end) = min(test_rep_temp, comm_rep_num);
                check_flag = true;
            end
            sel_counter = sel_counter+1;
%             rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
            rep_counter = rep_counter+comm_rep_num;
        end
        comm_selected = comm_all(I(1:sel_counter));
        rep_num_M = rep_num_M(1:sel_counter, :);
        test_rep_total = sum(rep_num_M(:, end));
        
        if check_flag == true
            % newborns of parents
            newborns_par(1:sl, 1:test_rep_total) = newborn_struct;
            rseed = randi(2^32-1, sl, test_rep_total, 'uint32');
            comm_type_counter = 0;
            rep_counter = 0;
            for i = 1:sel_counter
                comm_temp = pipette_SpikeTest_SP(comm_selected(i), newborn_struct, ...
                    const_struct, [spike_frac spike_test], rep_num_M(i, :), i);
                newborns(rep_counter+1:rep_counter+rep_num_M(i, 1)) ...
                    = comm_temp(1:rep_num_M(i, 1));
                rep_counter = rep_counter + rep_num_M(i, 1);
                for j = 1:sl
                    newborns_par(j, comm_type_counter+1:comm_type_counter+rep_num_M(i,end)) ...
                        = comm_temp(rep_num_M(i, 1)+1+(j-1)*rep_num_M(i,end) : rep_num_M(i, 1)+j*rep_num_M(i,end));
                end
                comm_type_counter = comm_type_counter+rep_num_M(i,end);
            end
            for i = 1:sl
                for j = 1:test_rep_total
                    newborns_par(i,j).rseed = rseed(i,j);
                end
            end
        end
    elseif abs(n - check_cycle + 1) < 0.1
%         check_flag = 1;           
        P_par = zeros(size(newborns_par));
        M0_par = zeros(size(newborns_par));
        H0_par = zeros(size(newborns_par));
        adults_par(1:sl, 1:test_rep_total) = comm_struct;
        parfor i = 1:test_rep_total*sl
            adults_par(i) = simu_one_comm(newborns_par(i), comm_struct, const_struct);
            P_par(i) = adults_par(i).P;
            M0_par(i) = adults_par(i).M_t(1);
            H0_par(i) = adults_par(i).H_t(1);
        end
        save([folder_name1 '/ParResults'],'P_par','M0_par','H0_par')
        newborns_off(1:test_rep_num*(sl+1), 1:test_rep_total) = newborn_struct;
        rseed = randi(2^32-1,test_rep_num*(sl+1), test_rep_total, 'uint32');
        sel_counter = 0;
        rep_counter = 0;
        for i = 1 : comm_type_num*comm_rep_num
            if rep_counter >= comm_type_num*comm_rep_num
                break
            end
            BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
            dil_factor = floor(BM/BM_target/(1-spike_frac));
            if dil_factor == 0
                continue
            end
            rep_num_temp = min([dil_factor-test_rep_num, comm_rep_num, N-rep_counter]);
            comm_temp = pipette_SpikeTest_SP(comm_all(I(i)), newborn_struct, ...
                const_struct, spike_frac, test_rep_num + rep_num_temp, i);
            newborns(rep_counter+1:rep_counter+rep_num_temp) ...
                = comm_temp(1:rep_num_temp);
            newborns_off(1 : test_rep_num, i) ...
                = transpose(comm_temp(rep_num_temp+1:rep_num_temp+test_rep_num));
            sel_counter = sel_counter+1;
            rep_counter = rep_counter+rep_num_temp;
        end
        comm_selected = comm_all(I(1:sel_counter));
        for i = sel_counter+1:test_rep_total
            newborns_off(1 : test_rep_num, i) ...
                = transpose(pipette_SpikeTest_SP(comm_all(I(i)), newborn_struct, ...
                const_struct, spike_frac, test_rep_num, i));
        end
        for i = 1:sl
            for j = 1:test_rep_total
                newborns_off(i*test_rep_num+1:(i+1)*test_rep_num, j) ...
                    = transpose(pipette_SpikeTest_SP(adults_par(i,j), newborn_struct, ...
                    const_struct, spike_test(i), test_rep_num, j));
            end
        end
        for i = 1:test_rep_num * (sl+1) * test_rep_total
            newborns_off(i).rseed = rseed(i);
        end
    elseif abs(n - check_cycle) < 0.1
        P_off = zeros(size(newborns_off));
        M0_off = zeros(size(newborns_off));
        H0_off = zeros(size(newborns_off));
        parfor i=1 : test_rep_total*test_rep_num*(sl+1)
            comm_temp = simu_one_comm(newborns_off(i), comm_struct, const_struct);
            P_off(i) = comm_temp.P;
            M0_off(i) = comm_temp.M_t(1);
            H0_off(i) = comm_temp.H_t(1);
        end
        check_flag = 0;
        save([folder_name1 '/OffResults'],'P_off','M0_off','H0_off')
    end
    if check_flag ==0
        rep_counter = 0;
        sel_counter = 0;
        %         comm_selected(1 : comm_type_num*comm_rep_num,1) = comm_struct;
        for i = 1 : comm_type_num*comm_rep_num
            if rep_counter >= comm_type_num*comm_rep_num
                break
            end
            dil_factor = floor((comm_all(I(i)).M_t(end)+comm_all(I(i)).H_t(end))/BM_target/(1-spike_frac));
            if dil_factor == 0
                continue
            end
            rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
            newborns(rep_counter+1:rep_counter+rep_num_temp) = pipette_SpikeTest_SP(comm_all(I(i)), newborn_struct, ...
                const_struct, spike_frac, rep_num_temp, i);
            sel_counter=sel_counter+1;
            rep_counter=rep_counter+rep_num_temp;
        end
        comm_selected=comm_all(I(1:sel_counter));
    end
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
    for ri = 1 : comm_type_num*comm_rep_num
        newborns(ri).rseed = rseed(ri);
    end
    
    % save the selected communities and the current state of the random number generator
    save([folder_name1 '/comm_selected'],'comm_selected');
    save([folder_name1 '/distrng'],'distrng');
    save([folder_name2 '/P_all'], 'P_all');
    save([folder_name2 '/Pn'],'Pn');
    n = n+1;
end
% save Newborns of the (C+1)th cycle
if ~exist('comm_all','dir')
    mkdir('comm_all')
end
save('comm_all/newborns', 'newborns')
save('comm_all/check_cycle','check_cycle')
