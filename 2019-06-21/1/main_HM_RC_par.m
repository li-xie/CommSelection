clear

C = 2e3; % total number of cycles
% reproducing method
% @pipette_RC: reproduce through pipetting, allowing both BM(0) and
%          phi_M(0) to fluctuate
% @cell_sort_RC: reproduce through cell-sorting, so that
%            the biomass in the Newborns are fixed BM(0)=BM_target and
%            phi_M(0)=phi_M(T) of the parent Adults
repro_method = @cell_sort_RC;

% the factor to be multiplied to R(0)
V = 1;
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
T0 = 17;
N = 100; % number of communities within a cycle
mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_type_num = 1; % minimal number of Adults allowed to reproduce
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4*V; % maximal number of cells in the community
t_bin = 0.05; % time step in the simulation
pcs=1e-15; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

% ancestral parameters shown in Table 1
gM_max = 0.7; % max growth rate of M
gH_max = 0.3; % max growth rate of H
dM = 3.5e-3; % death rate of M
dH = 1.5e-3; % death rate of H
fp_start = 0.1; % fp at the beginning of selection
r_B_start = 3;
c_BM_start = 1;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1/3;
K_HR = 1/5;
K_MB = 100;

% evolutionary bounds of the phenotypes
fp_Bound = 1; % fp is between 0 and 1
r_B_UBound = 6;
r_B_LBound = 1.5;
c_BM_UBound = 2;
c_BM_LBound = 0.5;


% R(0)=1 for each V
R0 = V;

% shuffle the random number generator state
rng('shuffle');
rseed=randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');

comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'r_B', zeros(max_popul,1), 'c_BM', zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target);

fp_init = zeros(max_popul, 1);
c_BM_init = zeros(max_popul, 1);
r_B_init = zeros(max_popul, 1);
M_L_init = zeros(max_popul, 1);
H_L_init=zeros(max_popul, 1);
% death_probability=zeros(max_popul,1);
M_counter = 60;
H_counter = BM_target - M_counter;
M_L_init(1 : M_counter) = 1;
H_L_init(1 : H_counter) = 1;
fp_init(1 : M_counter) = fp_start;
c_BM_init(1 : M_counter) = c_BM_start;
r_B_init(1 : H_counter) = r_B_start;


comm_all(1 : comm_type_num*comm_rep_num,1)...
    =struct('M_L', M_L_init,'H_L', H_L_init,'fp', fp_init,...
    'c_BM', c_BM_init, 'r_B', r_B_init,...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

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
        c_BM = comm_all(rep).c_BM;
        r_B = comm_all(rep).r_B;
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
            c_BM_M_L = c_BM(1:M_counter) .* M_L(1:M_counter);
            r_B_H_L = r_B(1:H_counter) .* H_L(1:H_counter);
            para = [gM_max; gH_max; c_RM; c_RH; K_MR; K_HR; K_MB; sum(M_L); sum(H_L); sum(c_BM_M_L); sum(r_B_H_L)];
            fhandle=@(t,y) chem_conc_jacobian_RC(t, y , para);
            options=odeset('Jacobian',fhandle,'RelTol',1e-5);
            [tx,y]=ode23s(@(t,y) chem_conc_RC(t, y, para), [0 t_bin], [B(dt-1); R(dt-1)], options);
            if ~isreal(y)
                error('imaginary value')
            end
            % matrice where the row number is the number of cells and colomn
            % number is the B(t) and R(t) within the time step
            BN = y(:,1)/K_MB;
            RN_M = y(:,2)/K_MR;
            RN_H = y(:,2)/K_HR;
            M_coef = BN./(RN_M+BN).*RN_M./(RN_M+1) ...
                + RN_M./(RN_M+BN).*BN./(BN+1);
            gMdt = trapz(tx, M_coef) * gM_max.* (1 - fp(1:M_counter));
            H_coef = RN_H ./(RN_H + 1);
            gHdt = trapz(tx, H_coef) * gH_max;
            
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
            M_t(dt) = sum(M_L);
            
            death_probability = rand(max_popul,1);
            H_L( death_probability < dH * t_bin) = 0;
            H_t(dt) = sum(H_L);
            
            div_idx = find(M_L >= 2);
            div_length = length(div_idx);
            if div_length > 0
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                fp(M_counter+1 : M_counter+div_length)...
                    = fp(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                c_BM(M_counter+1 : M_counter+div_length)...
                    = c_BM(div_idx) .*mut_multiplier;
                
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                fp(div_idx) = fp(div_idx) .* mut_multiplier;
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                c_BM(div_idx)=c_BM(div_idx) ./mut_multiplier;
                
                M_L(M_counter+1:M_counter+div_length) = M_L(div_idx)/2;
                M_L(div_idx) = M_L(div_idx)/2;
                M_counter = M_counter + div_length;
                
                % if the phenotypes exceed the bounds, cap them at the bounds
                fp(fp > fp_Bound) = fp_Bound;
                c_BM(c_BM < c_BM_LBound) = c_BM_LBound;
                c_BM(c_BM > c_BM_UBound) = c_BM_UBound;
            end
            
            div_idx = find(H_L >= 2);
            div_length = length(div_idx);
            if div_length > 0
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                r_B(H_counter+1 : H_counter+div_length)...
                    =r_B(div_idx) .* mut_multiplier;
                
                mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
                r_B(div_idx) = r_B(div_idx) .* mut_multiplier;
                
                H_L(H_counter+1 : H_counter+div_length) = H_L(div_idx)/2;
                H_L(div_idx) = H_L(div_idx)/2;
                H_counter = H_counter + div_length;
                
                r_B(r_B < r_B_LBound) = r_B_LBound;
                r_B(r_B > r_B_UBound) = r_B_UBound;
            end
            
        end
        % find the indice of the cells that are viable (positive growth 
        % rates and K smaller than the singular value)
        temp1 = find((M_L > pcs));
        temp2 = find((H_L > pcs));
        M_t(end) = sum(M_L(temp1));
        H_t(end) = sum(H_L(temp2));
        M_counter = length(temp1);
        comm_rep.M_L(1:M_counter) = M_L(temp1);
        H_counter = length(temp2);
        comm_rep.H_L(1:H_counter) = H_L(temp2);
        
        comm_rep.fp(1:M_counter) = fp(temp1);
        comm_rep.c_BM(1:M_counter) = c_BM(temp1);
        comm_rep.r_B(1:H_counter) = r_B(temp2);
        
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
    save([folder_name2 '/P_all'], 'P_all');
end

% if ~exist('comm_all','dir')
%     mkdir('comm_all')
% end
% for i=1:comm_type_num
%     comm=comm_all((i-1)*comm_rep_num+1:i*comm_rep_num);
%     save(['comm_all/comm' num2str(i*comm_rep_num)],'comm');
% end
