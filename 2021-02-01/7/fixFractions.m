% reproduce Adult communities into Newborns with fixed phi_M(0) = phi_M(T)
% of the parent Adult community
% comm_select: Adults to be reproduced
% const_struct: a structure to pass constants
% dil_factor: dilution factor
% rep_counter: current number of Newborn communities.
% parentnum: the rank of the Adult community
function comm_rep = fixFractions(comm_selected,comm_struct, const_struct, dil_factor, rep_counter, parentnum)
% maximal number of offspring community from one Adult
comm_rep_num = const_struct.comm_rep_num;
% minimal number of Adults allowed to reproduce
comm_type_num = const_struct.comm_type_num;
% comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
%     'M_t',zeros(t_binnum,1),'H_t',zeros(t_binnum,1),'R',zeros(t_binnum,1),'B',zeros(t_binnum,1),...
%     'P',0,'parentnum',0,'rseed',uint32(0));comm_rep(1:comm_rep_num,1)=comm_struct;
comm_rep(1:comm_rep_num,1)=comm_struct;
M_counter = nnz(comm_selected.M_L);
H_counter = nnz(comm_selected.H_L);
M_L = comm_selected.M_L(1:M_counter);
H_L = comm_selected.H_L(1:H_counter);
fp = comm_selected.fp(1:M_counter);

% calculate the H:M biomass in the parent Adult
HMRatio=sum(H_L)/sum(M_L);
% randomly sort M cells into dil_factor Newborns
rand_temp1 = ceil(rand(M_counter,1) * dil_factor);

% randomly permute H cells
rand_idx2=randperm(H_counter);
H_L_rand=H_L(rand_idx2);
% initialize the partition indice for H cells
par_idx2=[0; zeros(dil_factor,1)];
% if the Adult can generate more Newborns than comm_rep_num, keep only
% comm_rep_num
if dil_factor >= comm_rep_num
    for i = 1 : comm_rep_num
        % sort all the M cells with random number i into the ith Newborn
        temp_idx1 = find(rand_temp1 == i);
        M_num = length(temp_idx1);
        if M_num >= 1
            comm_rep(i).M_L(1 : M_num) = M_L(temp_idx1);
            comm_rep(i).fp(1 : M_num) = fp(temp_idx1);
        end
        % assign H cells into each Newborn according to the M biomass in
        % the Newborn so that the H:M is the same as that of its parent
        % Adult community
        % the H cells that haven't been assigned to Newborns
        H_L_remain = H_L_rand(par_idx2(i)+1:H_counter);
        % calculate the accumulative biomass
        H_remain_accu = cumsum(H_L_remain);
        % find the group of H cells so that their biomass is the closest to without exceeding 
        % M_biomass*HMRatio 
        [~, idx_temp] = min( abs(H_remain_accu - sum(M_L(temp_idx1))*HMRatio) );
        if H_remain_accu(idx_temp)-sum(M_L(temp_idx1))*HMRatio > 0
            par_idx2(i+1) = idx_temp - 1 + par_idx2(i);
        else
            par_idx2(i+1) = idx_temp + par_idx2(i);
        end
        H_num = par_idx2(i+1) - par_idx2(i);
        if H_num >= 1
            comm_rep(i).H_L(1 : H_num) = H_L_rand(par_idx2(i)+1 : par_idx2(i+1));
        end
        comm_rep(i).parentnum = parentnum;
    end
    if rep_counter+comm_rep_num > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    end
else
    for i= 1 : dil_factor
        % sort all the M cells with random number i into the ith Newborn
        temp_idx1 = find(rand_temp1 == i);
        M_num = length(temp_idx1);
        if M_num >= 1
            comm_rep(i).M_L(1 : M_num) = M_L(temp_idx1);
            comm_rep(i).fp(1 : M_num) = fp(temp_idx1);
        end
        % assign H cells into each Newborn according to the M biomass in
        % the Newborn so that the H:M is the same as that of its parent
        % Adult community
        % the H cells that haven't been assigned to Newborns
        H_L_remain = H_L_rand(par_idx2(i)+1:H_counter);
        % calculate the accumulative biomass
        H_remain_accu = cumsum(H_L_remain);
        % find the group of H cells so that their biomass is the closest to without exceeding 
        % M_biomass*HMRatio
        [~, idx_temp] = min( abs(H_remain_accu-sum(M_L(temp_idx1))*HMRatio) );
        if H_remain_accu(idx_temp)-sum(M_L(temp_idx1))*HMRatio > 0
            par_idx2(i+1) = idx_temp-1 + par_idx2(i);
        else
            par_idx2(i+1) = idx_temp + par_idx2(i);
        end
        H_num = par_idx2(i+1) - par_idx2(i);
        if H_num >= 1
            comm_rep(i).H_L(1:H_num) = H_L_rand(par_idx2(i)+1 : par_idx2(i+1));
        end
        comm_rep(i).parentnum = parentnum;
    end
    if rep_counter+dil_factor > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    else
        comm_rep = comm_rep(1 : dil_factor);
    end
end