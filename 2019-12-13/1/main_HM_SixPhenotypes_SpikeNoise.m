clear

C = 2; % total number of cycles
% upper bound for gH_max. for gH_max_Bound = 0.3, set V = 1.
gH_max_Bound = 0.8;
% the factor for amount of R(0)
V = 10;
% minimal number of Adults allowed to reproduce. comm_type_num = 1 for the
% top-dog strategy, comm_type_num = n for the top n% strategy.
comm_type_num = 1; 
% mutation rate corresponding to effective mutation rate of 2e-3
% to turn off the mutation, set mut_rate=0.
mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
% reproducing method
% @pipette_SP: reproduce through pipetting, allowing both BM(0) and
%          phi_M(0) to fluctuate
% @cell_sort_SP: reproduce through cell-sorting, so that
%            the biomass in the Newborns are fixed BM(0)=BM_target and
%            phi_M(0)=phi_M(T) of the parent Adults
repro_method = @cell_sort_SP;

N = 100; % number of communities within a cycle
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
T0 = 17;
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4*V; % maximal number of cells in the community
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

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target);

fp_init = zeros(max_popul, 1);
gM_max_init = zeros(max_popul, 1);
K_MB_init = zeros(max_popul, 1);
K_MR_init = zeros(max_popul, 1);
gH_max_init = zeros(max_popul, 1);
M_L_init = zeros(max_popul, 1);
K_HR_init = zeros(max_popul, 1);
H_L_init=zeros(max_popul, 1);
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


comm_all(1 : comm_type_num*comm_rep_num,1)...
    =struct('M_L', M_L_init,'H_L', H_L_init,'fp', fp_init,...
    'gM_max', gM_max_init, 'K_MB', K_MB_init, 'K_MR', K_MR_init,...
    'gH_max', gH_max_init, 'K_HR', K_HR_init,...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
rseed=randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
for ri=1:comm_type_num*comm_rep_num
    comm_all(ri).rseed=rseed(ri);
end


for n = 1 : C
    % create a folder Cn to save the results of the nth cycle
    folder_name1 = ['C' num2str(n)];
    if ~exist(folder_name1, 'dir')
        mkdir(folder_name1)
    end
    folder_name2=['C' num2str(n) '/comm_all'];
    if ~exist(folder_name2, 'dir')
        mkdir(folder_name2)
    end
    newborns = newborn_trim(comm_all, BM_target, pcs);
    save([folder_name2 '/newborns'],'newborns');
    % rep is the index of communities within one cycle
    parfor rep = 1:comm_type_num*comm_rep_num
        rnodeseed = comm_all(rep).rseed;
        rng(rnodeseed, 'twister');
        comm_rep = comm_struct;
        % copy the data from the structure with Newborn's configuration
        M_L = comm_all(rep).M_L;
        H_L = comm_all(rep).H_L;
        fp = comm_all(rep).fp;
        gM_max = comm_all(rep).gM_max;
        K_MB = comm_all(rep).K_MB;
        K_MR = comm_all(rep).K_MR;
        gH_max = comm_all(rep).gH_max;
        K_HR = comm_all(rep).K_HR;
        M_t = zeros(t_binnum+1, 1);
        H_t = zeros(t_binnum+1, 1);
        %  temporary column vectors for M and H's biomass
        M_LTemp = zeros(max_popul,1);
        H_LTemp = zeros(max_popul,1);
        % column vectors to store B and R's concentrations at each time
        % step
        B = zeros(t_binnum+1, 1);
        R = zeros(t_binnum+1, 1);
        P = 0;
        %         p0=sum(fp.*M_L./(1-fp));
        M_counter = nnz(M_L>pcs);
        H_counter = nnz(H_L>pcs);
        M_t(1) = sum(M_L);
        H_t(1) = sum(H_L);
        R(1) = R0;
        % perform simulation through t_binnum time steps
        for dt = 2 : t_binnum+1
            gM_maxM_L = gM_max(1:M_counter) .* M_L(1:M_counter);
            gH_maxH_L = gH_max(1:H_counter) .* H_L(1:H_counter);
            paras=struct('gM_maxM_L',gM_maxM_L,'gH_maxH_L',gH_maxH_L,...
                'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH,...
                'K_MB', K_MB(1:M_counter),'K_MR',K_MR(1:M_counter),...
                'K_HR',K_HR(1:H_counter));
            fhandle=@(t,y) chem_conc_jacobian_SP(t, y , paras);
            options=odeset('Jacobian',fhandle,'RelTol',1e-5);
            [tx,y]=ode23s(@(t,y) chem_conc_SP(t, y, paras), [0 t_bin], [B(dt-1); R(dt-1)], options);
            if ~isreal(y)
                error('imaginary value')
            end
            % matrice where the row number is the number of cells and colomn
            % number is the B(t) and R(t) within the time step
            BN_mat = (1./K_MB(1:M_counter)) * y(:,1)';
            RN_M_mat = (1./K_MR(1:M_counter)) * y(:,2)';
            M_coef = BN_mat./(RN_M_mat+BN_mat).*RN_M_mat./(RN_M_mat+1) ...
                + RN_M_mat./(RN_M_mat+BN_mat).*BN_mat./(BN_mat+1);
            gMdt = trapz(tx, M_coef, 2) .* gM_max(1:M_counter) .* (1 - fp(1:M_counter));
            RN_H_mat = (1./K_HR(1:H_counter)) * y(:,2)';
            H_coef = RN_H_mat ./(RN_H_mat+1);
            gHdt = trapz(tx, H_coef, 2) .* gH_max(1:H_counter);
            
            M_LTemp(1:M_counter) = exp(gMdt) .* M_L((1:M_counter));
            H_LTemp(1:H_counter)=exp(gHdt) .* H_L(1:H_counter);
            P = sum(fp(1:M_counter) .* (M_LTemp(1:M_counter) - M_L(1:M_counter))...
                ./(1-fp(1:M_counter)) ) + P;
            M_L(1:M_counter) = M_LTemp(1:M_counter);
            H_L(1:H_counter) = H_LTemp(1:H_counter);
            B(dt) = y(end, 1);
            R(dt) = y(end, 2);
            
            death_probability = rand(max_popul,1);
            M_L(death_probability < dM * t_bin) = 0;
            fp( death_probability < dM * t_bin) = 0;
            gM_max( death_probability < dM * t_bin) = 0;
            M_t(dt) = sum(M_L);
            
            death_probability = rand(max_popul,1);
            H_L( death_probability < dH * t_bin) = 0;
            gH_max( death_probability < dH * t_bin) = 0;
            H_t(dt) = sum(H_L);
            
            div_idx = find(M_L >= 2);
            div_length = length(div_idx);
            if div_length > 0
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                fp(M_counter+1 : M_counter+div_length)...
                    = fp(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_MB(M_counter+1 : M_counter+div_length)...
                    = K_MB(div_idx) ./mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_MR(M_counter+1 : M_counter+div_length)...
                    = K_MR(div_idx) ./mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                gM_max(M_counter+1 : M_counter+div_length)...
                    =gM_max(div_idx) .* mut_multiplier;
                
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                fp(div_idx) = fp(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_MB(div_idx)=K_MB(div_idx) ./mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_MR(div_idx)=K_MR(div_idx) ./mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                gM_max(div_idx)=gM_max(div_idx) .* mut_multiplier;
                
                M_L(M_counter+1:M_counter+div_length) = M_L(div_idx)/2;
                M_L(div_idx) = M_L(div_idx)/2;
                M_counter = M_counter + div_length;
                
                % if the phenotypes exceed the bounds, cap them at the bounds
                fp(fp > fp_Bound) = fp_Bound;
                K_MB((K_MB < K_MB_Bound)&(K_MB > pcs)) = K_MB_Bound;
                K_MR((K_MR < K_MR_Bound)&(K_MR > pcs)) = K_MR_Bound;
                gM_max(gM_max > gM_max_Bound) = gM_max_Bound;
            end
            
            div_idx = find(H_L >= 2);
            div_length = length(div_idx);
            if div_length > 0
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                gH_max(H_counter+1 : H_counter+div_length)...
                    =gH_max(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_HR(H_counter+1 : H_counter+div_length)...
                    = K_HR(div_idx) ./mut_multiplier;
                
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                gH_max(div_idx) = gH_max(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                mut_multiplier = max(mut_multiplier, pcs);
                K_HR(div_idx) = K_HR(div_idx) ./mut_multiplier;
                
                H_L(H_counter+1 : H_counter+div_length) = H_L(div_idx)/2;
                H_L(div_idx) = H_L(div_idx)/2;
                H_counter = H_counter + div_length;
                
                gH_max(gH_max > gH_max_Bound) = gH_max_Bound;
                K_HR((K_HR < K_HR_Bound)&(K_HR > pcs)) = K_HR_Bound;
            end
            
        end
        % find the indice of the cells that are viable (positive growth 
        % rates and K smaller than the singular value)
        temp1 = find((gM_max > pcs) & (K_MB < K_singular) & (K_MR < K_singular));
        temp2 = find((gH_max > pcs) & (K_HR < K_singular));
        M_t(end) = sum(M_L(temp1));
        H_t(end) = sum(H_L(temp2));
        M_counter = length(temp1);
        comm_rep.M_L(1:M_counter) = M_L(temp1);
        H_counter = length(temp2);
        comm_rep.H_L(1:H_counter) = H_L(temp2);
        
        comm_rep.fp(1:M_counter) = fp(temp1);
        comm_rep.gM_max(1:M_counter) = gM_max(temp1);
        comm_rep.K_MB(1:M_counter) = K_MB(temp1);
        comm_rep.K_MR(1:M_counter) = K_MR(temp1);
        comm_rep.gH_max(1:H_counter) = gH_max(temp2);
        comm_rep.K_HR(1:H_counter) = K_HR(temp2);
        
        comm_rep.M_t = M_t;
        comm_rep.H_t = H_t;
        comm_rep.B = B;
        comm_rep.R = R;
        comm_rep.P = P;
        comm_rep.parentnum = comm_all(rep).parentnum;
        comm_rep.rseed = comm_all(rep).rseed;
        
        comm_all(rep)=comm_rep;
    end
    distrng=rng;
    P_all = [comm_all.P];
    [~, I] = sort([comm_all.P], 'descend');
    % I=randperm(comm_type_num*comm_rep_num);
    comm_all_sorted = comm_all(I);
    rep_counter = 0;
    sel_counter = 0;
    comm_selected(1 : comm_type_num*comm_rep_num,1) = comm_struct;
    for i = 1 : comm_type_num*comm_rep_num
        if rep_counter >= comm_type_num*comm_rep_num
            break
        end
        dil_factor = floor((comm_all_sorted(i).M_t(t_binnum+1)+comm_all_sorted(i).H_t(t_binnum+1))/BM_target);
        if dil_factor == 0
            continue
        end
        rep_num_temp = min(dil_factor, comm_rep_num);
        comm_all_idx = min(comm_type_num*comm_rep_num, rep_counter+rep_num_temp);
        comm_all(rep_counter+1:comm_all_idx) = repro_method(comm_all_sorted(i), comm_struct, const_struct, dil_factor, rep_counter, i);
        sel_counter=sel_counter+1;
        comm_selected(sel_counter)=comm_all_sorted(i);
        rep_counter=rep_counter+rep_num_temp;
    end
    comm_selected=comm_selected(1:sel_counter);
    
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
    for ri = 1 : comm_type_num*comm_rep_num
        comm_all(ri).rseed = rseed(ri);
    end
    % save the selected communities and the current state of the random number generator 
    save([folder_name1 '/comm_selected'],'comm_selected');
    save([folder_name1 '/distrng'],'distrng');
%     save([folder_name2 '/P_all'], 'P_all');
end

