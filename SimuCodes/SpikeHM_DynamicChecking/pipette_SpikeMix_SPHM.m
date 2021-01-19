% reproduce Adult communities to estimate heritability under different substitution fraction,
% when six phenotypes could mutate. Modified to use multinomial distribution to draw 
% number of cells into each Newborn. Additionally, randomly pick one H cell.

% comm_selected: Adult to be reproduced. 1x1 structure.
% newborn_struct: 1x1 structure to initialize Newborns
% const_struct: a structure to pass constants

% spike_frac_M: 1 by sl+1 row vector that specifies the fraction of Newborns' biomass to
% be substituted with H biomass isolated from the same lineage from the previous cycle.
% spike_frac_M(2:end) have the same number

% rep_num_M: 1 by sl+1 vector that specifies the number of Newborns corresponding to each
% substitution fraction specified in spike_frac_M.

% H_isolate_in: 1x1 structure of the H cell used for substitution
% parentnum: the rank of the Adult community. Used to keep track of lineages

%comm_rep is the Newborns reproduced from comm_selected
% H_isolate_out is the H cell randomly picked from the Adult after reproduction
function [comm_rep, H_isolate_out, M_isolate_out] = pipette_SpikeMix_SPHM(comm_selected, newborn_struct, ...
    const_struct, spike_frac_M, rep_num_M, H_isolate_in, M_isolate_in, spike_clone_num, parentnum)
BM_target = const_struct.BM_target;
pcs = const_struct.pcs;
% comm_rep is a structure row array
comm_rep(1,1:sum(rep_num_M)) = newborn_struct;
gH_max_spike = H_isolate_in.gH_max;
K_HR_spike = H_isolate_in.K_HR;
H_L_spike = H_isolate_in.H_L;
gM_max_spike = M_isolate_in.gM_max;
K_MR_spike = M_isolate_in.K_MR;
K_MB_spike = M_isolate_in.K_MB;
fp_spike = M_isolate_in.fp;
M_L_spike = M_isolate_in.M_L;
if length(M_L_spike) ~= spike_clone_num...
        || length(H_L_spike) ~= spike_clone_num
    save('debug_data')
    error ('number of HM isolate is not right')
end
temp_idx = find(M_L_spike > pcs);
gM_max_spike = gM_max_spike(temp_idx);
K_MR_spike = K_MR_spike(temp_idx);
K_MB_spike = K_MB_spike(temp_idx);
fp_spike = fp_spike(temp_idx);
M_L_spike = M_L_spike(temp_idx);
M_L_mean = mean(M_L_spike);
M_clone_num = length(M_L_spike);

temp_idx = find(H_L_spike > pcs);
gH_max_spike = gH_max_spike(temp_idx);
K_HR_spike = K_HR_spike(temp_idx);
H_L_spike = H_L_spike(temp_idx);
H_L_mean = mean(H_L_spike);
H_clone_num = length(H_L_spike);

M_counter = nnz(comm_selected.M_L);
H_counter = nnz(comm_selected.H_L);
M_L = comm_selected.M_L(1:M_counter);
H_L = comm_selected.H_L(1:H_counter);
fp = comm_selected.fp(1:M_counter);
gM_max = comm_selected.gM_max(1:M_counter);
gH_max = comm_selected.gH_max(1:H_counter);
K_MB = comm_selected.K_MB(1:M_counter);
K_MR = comm_selected.K_MR(1:M_counter);
K_HR = comm_selected.K_HR(1:H_counter);

all_L = [M_L; H_L];
all_fp = [fp; -1*ones(H_counter,1)];
all_gMax = [gM_max; gH_max];
all_KR = [K_MR; K_HR];
all_KB = [K_MB; -1*ones(H_counter,1)];
BM = sum(all_L);
% randomize all cells before distributing into Newborns
rand_idx = randperm(M_counter + H_counter);
all_L_rand = all_L(rand_idx);
all_fp_rand = all_fp(rand_idx);
all_gMax_rand = all_gMax(rand_idx);
all_KR_rand = all_KR(rand_idx);
all_KB_rand = all_KB(rand_idx);

if ~isempty(spike_frac_M)
    % calculate the number of cells into each Newborn from multinomial distribution. 
    % mn_p is the probability of drawing a cell for each Newborn. e.g., if 30% of the Newborn
    % will be substituted with H, the probability that a cell goes to this Newborn is
    % (1-30%)*BM_target/BM, where BM is the total biomass of the Adult
    % p_temp has rep_num_M(end) rows and sl columns 
    HM_ind = zeros(sum(rep_num_M), 1);
    p_temp = ones(rep_num_M(end), 1) * (BM_target*(1-abs(spike_frac_M(2:end)))/BM);
    mn_p = [BM_target*(1-abs(spike_frac_M(1)))/BM*ones(1,rep_num_M(1)), transpose(p_temp(:))];
    if sum(mn_p)>1 || nnz(mn_p<0)>0
        save('debug_data')
        error ('mistake in number of Newborns')
    end
    mn_p = [mn_p, 1-sum(mn_p)];
    % draw the number of cells to each Newborn from a multinomial distribution
    num_cells = mnrnd(M_counter + H_counter, mn_p);
    num_cells = num_cells(1:end-1);
    % H_spike is the number of H cells assigned to each Newborn. Each row corresponds to
    % number of H cells spiked into each Newborn
    HM_spike = zeros(sum(rep_num_M), 1);
    temp_counter = 0;
    for s = 1:length(spike_frac_M)
        if spike_frac_M(s)>0
            HM_ind(temp_counter+1:temp_counter+rep_num_M(s)) = 1;
            % draw the number of H cells assigned to each Newborn from a Poisson distribution
            HM_spike(temp_counter+1:temp_counter+rep_num_M(s)) ...
                = poissrnd(BM_target*spike_frac_M(s)/H_L_mean, [rep_num_M(s), 1]);
        elseif spike_frac_M(s)<0
            HM_ind(temp_counter+1:temp_counter+rep_num_M(s)) = -1;
            HM_spike(temp_counter+1:temp_counter+rep_num_M(s)) ...
            = poissrnd(-BM_target*spike_frac_M(s)/M_L_mean, [rep_num_M(s), 1]);
        elseif abs(spike_frac_M(s))<pcs
            HM_ind(temp_counter+1:temp_counter+rep_num_M(s)) = 0;
            HM_spike(temp_counter+1:temp_counter+rep_num_M(s)) = 0;
        else
            save('debug_data')
            error('error in assigning spiking numbers')
        end            
        temp_counter = temp_counter + rep_num_M(s);
    end
end
if sum(num_cells) > M_counter+H_counter
    error('not enough cells')
else
    H_isolate_out.H_L = -ones(spike_clone_num, 1);
    H_isolate_out.gH_max = -ones(spike_clone_num, 1);
    H_isolate_out.K_HR = -ones(spike_clone_num, 1);
    M_isolate_out.M_L = -ones(spike_clone_num, 1);
    M_isolate_out.gM_max = -ones(spike_clone_num, 1);
    M_isolate_out.K_MR = -ones(spike_clone_num, 1);
    M_isolate_out.fp = -ones(spike_clone_num, 1);
    M_isolate_out.K_MB = -ones(spike_clone_num, 1);
    % randomly pick a H cell for substitution during the next cycle
    remain_start = sum(num_cells)+1;
    H_idx = find(all_fp_rand(remain_start:end)<-pcs, spike_clone_num);
    if ~isempty(H_idx)
        t_l = length(H_idx);
        H_idx = H_idx+remain_start-1;
        H_isolate_out.H_L(1:t_l) = all_L_rand(H_idx);
        H_isolate_out.gH_max(1:t_l) = all_gMax_rand(H_idx);
        H_isolate_out.K_HR(1:t_l) = all_KR_rand(H_idx);
        if nnz(H_isolate_out.gH_max(1:t_l) > const_struct.gH_max_Bound+pcs)>0
            save('debug_data')
            error('M isolate is mistaken for H isolate')
        end
    else
        H_isolate_out = H_isolate_in;
    end
    M_idx = find(all_fp_rand(remain_start:end)>-pcs, spike_clone_num);
    if ~isempty(M_idx)
        t_l = length(M_idx);
        M_idx = M_idx+remain_start-1;
        M_isolate_out.M_L(1:t_l) = all_L_rand(M_idx);
        M_isolate_out.gM_max(1:t_l) = all_gMax_rand(M_idx);
        M_isolate_out.K_MR(1:t_l) = all_KR_rand(M_idx);
        M_isolate_out.fp(1:t_l) = all_fp_rand(M_idx);
        M_isolate_out.K_MB(1:t_l) = all_KB_rand(M_idx);
        if nnz(M_isolate_out.K_MB(1:t_l) < 0)>0
            save('debug_data')
            error('H isolate is mistaken for M isolate')
        end
    else
        M_isolate_out = M_isolate_in;
    end
end
temp_counter = 0;
for i = 1 : sum(rep_num_M)
    L_temp = all_L_rand(temp_counter+1:temp_counter+num_cells(i));
    fp_temp = all_fp_rand(temp_counter+1:temp_counter+num_cells(i));
    gMax_temp = all_gMax_rand(temp_counter+1:temp_counter+num_cells(i));
    KR_temp = all_KR_rand(temp_counter+1:temp_counter+num_cells(i));
    KB_temp = all_KB_rand(temp_counter+1:temp_counter+num_cells(i));
    M_num = nnz(fp_temp>-0.1);
    H_num = nnz(fp_temp<-0.1);
    temp_counter = temp_counter + num_cells(i);
    if M_num >= 1
        comm_rep(i).M_L(1:M_num) = L_temp(fp_temp>-0.1);
        comm_rep(i).fp(1:M_num) = fp_temp(fp_temp>-0.1);
        comm_rep(i).gM_max(1:M_num) = gMax_temp(fp_temp>-0.1);
        comm_rep(i).K_MB(1:M_num) = KB_temp(fp_temp>-0.1);
        comm_rep(i).K_MR(1:M_num) = KR_temp(fp_temp>-0.1);
        if ~isempty(find(((comm_rep(i).K_MR < (const_struct.K_MR_Bound-pcs))...
                & (comm_rep(i).K_MR >0)), 1))
            save('debug_data')
            error('H cell is mistaken for M cell')
        end
    end
    if H_num >= 1
        comm_rep(i).H_L(1:H_num) = L_temp(fp_temp<-0.1);
        comm_rep(i).gH_max(1:H_num) = gMax_temp(fp_temp<-0.1);
        comm_rep(i).K_HR(1:H_num) = KR_temp(fp_temp<-0.1);
        if ~isempty(find(comm_rep(i).gH_max > (const_struct.gH_max_Bound+pcs), 1))
            save('debug_data')
            error('M cell is mistaken for H cell')
        end
    end
    % supplement with H or M cells with identical phenotype from the previous cycle of the same
    % lineage
    if HM_ind(i) > pcs
        if HM_spike(i) >= 1
            spike_idx = randsample(H_clone_num, HM_spike(i), true);
            comm_rep(i).H_L(1+H_num : HM_spike(i)+H_num) = H_L_spike(spike_idx);
            comm_rep(i).gH_max(1+H_num : HM_spike(i)+H_num) = gH_max_spike(spike_idx);
            comm_rep(i).K_HR(1+H_num : HM_spike(i)+H_num) = K_HR_spike(spike_idx);
        end
    elseif HM_ind(i) < -pcs
        if HM_spike(i) >= 1
            spike_idx = randsample(M_clone_num, HM_spike(i), true);
            comm_rep(i).M_L(1+M_num : HM_spike(i)+M_num) = M_L_spike(spike_idx);
            comm_rep(i).gM_max(1+M_num : HM_spike(i)+M_num) = gM_max_spike(spike_idx);
            comm_rep(i).K_MR(1+M_num : HM_spike(i)+M_num) = K_MR_spike(spike_idx);
            comm_rep(i).K_MB(1+M_num : HM_spike(i)+M_num) = K_MB_spike(spike_idx);
            comm_rep(i).fp(1+M_num : HM_spike(i)+M_num) = fp_spike(spike_idx);
        end
    elseif abs(HM_ind(i)) < pcs
        if abs(HM_spike(i)) >= 1
            save('debug_data')
            error('spiking numbers do not match spiking type')
        end
    end
    comm_rep(i).parentnum = parentnum;
end




