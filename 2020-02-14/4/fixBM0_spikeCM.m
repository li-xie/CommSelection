% reproduce Adult communities into Newborns with fixed BM(0)=BM_target
% comm_select: Adults to be reproduced
% const_struct: a structure to pass constants
% dil_factor: dilution factor
% rep_counter: current number of Newborn communities.
% parentnum: the rank of the Adult community
function comm_rep = fixBM0_spikeCM(comm_selected, comm_struct, const_struct,...
    dil_factor, spike_frac, rep_counter, parentnum)
% maximal number of offspring community from one Adult
comm_rep_num = const_struct.comm_rep_num;
% minimal number of Adults allowed to reproduce
comm_type_num = const_struct.comm_type_num;
BM_target = const_struct.BM_target;
M_Biomass_spike = BM_target*spike_frac;
% comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
%     'M_t',zeros(t_binnum,1),'H_t',zeros(t_binnum,1),'R',zeros(t_binnum,1),'B',zeros(t_binnum,1),...
%     'P',0,'parentnum',0,'rseed',uint32(0));
% initialize for the Newborn communities
comm_rep(1 : comm_rep_num,1) = comm_struct;

M_counter = nnz(comm_selected.M_L);
H_counter = nnz(comm_selected.H_L);
M_L = comm_selected.M_L(1:M_counter);
H_L = comm_selected.H_L(1:H_counter);
fp = comm_selected.fp(1:M_counter);
% concatenate all the cells together
all_L = [M_L; H_L];
% if fp==2, it's a H cell
all_fp = [fp; 2*ones(H_counter,1)];
% randomly permute cells
rand_idx = randperm(M_counter + H_counter);
all_L_rand = all_L(rand_idx);
all_fp_rand = all_fp(rand_idx);

% distribute cells randomly into a Newborn until the biomass reaches
% BM_target without exceeding BM_target
all_accu = cumsum(all_L_rand);
% partition index for different Newborns 
par_idx = [0; zeros(dil_factor,1)];
% Biomass distributed from Adult communities
BM_dist = BM_target - M_Biomass_spike;
for i = 1 : dil_factor
    [~, par_idx(i+1)]=min(abs(all_accu - i*BM_dist));
end
% temporary matrix for biomass of M and H cells
L_temp=zeros(BM_dist,1);
% temporary matrix for fp of M and H cells. if fp==2, it's a H cell
fp_temp=zeros(BM_dist,1);
% if the Adult can generate more Newborns than comm_rep_num, keep only
% comm_rep_num
if dil_factor >= comm_rep_num
    remain_start = par_idx(comm_rep_num+1)+1;
    M_idx = find(all_fp_rand(remain_start:end) <= 1) + remain_start - 1;
    M_L_remain = all_L_rand(M_idx);
    fp_remain = all_fp_rand(M_idx);
    if nnz(fp_remain>1)>0
        save('debug_data')
        error('H is mistaken for M')
    end
    M_accu = cumsum(M_L_remain);
    % partition index for different Newborns
    M_par_idx = [0; zeros(comm_rep_num,1)];
    % Biomass distributed from Adult communities
    for i = 1 : comm_rep_num
        [~, M_par_idx(i+1)] = min(abs(M_accu - i*M_Biomass_spike));
    end
    for i = 1:comm_rep_num
        % collect all the biomass and fp values between two partition indice
        L_temp(1:par_idx(i+1)-par_idx(i)) = all_L_rand(par_idx(i)+1:par_idx(i+1));
        fp_temp(1:par_idx(i+1)-par_idx(i)) = all_fp_rand(par_idx(i)+1:par_idx(i+1));
        M_idx = find(fp_temp(1:par_idx(i+1)-par_idx(i)) <= 1);
        H_idx = find(fp_temp(1:par_idx(i+1)-par_idx(i)) > 1);
        M_num = length(M_idx);
        M_spike_num = M_par_idx(i+1) - M_par_idx(i);
        H_num = length(H_idx);
        if M_num >= 1
            comm_rep(i).M_L(1:M_num) = L_temp(M_idx);
            comm_rep(i).fp(1:M_num) = fp_temp(M_idx);
            comm_rep(i).M_L(1+M_num : M_spike_num+M_num) ...
                = M_L_remain(M_par_idx(i)+1 : M_par_idx(i+1));
            comm_rep(i).fp(1+M_num : M_spike_num+M_num) ...
                = fp_remain(M_par_idx(i)+1 : M_par_idx(i+1));
        end
        if H_num >= 1
            comm_rep(i).H_L(1:H_num) = L_temp(H_idx);
%             comm_rep(i).H_L(1+H_num : H_spike+H_num) = H_L(1);
        end
        comm_rep(i).parentnum = parentnum;
    end
    if rep_counter+comm_rep_num > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1:comm_type_num*comm_rep_num-rep_counter);
    end
    
% if the Adult generates less Newborns than comm_rep_num, keep them all
else
    error('not enough cells for comm_rep_num Newborns')
end
    
    
