% reproduce Adult communities for spiking test, when six phenotypes could
% mutate. Modified to use multinomial distribution to draw number of cells 
% randomly select one H cell for future spiking 

% comm_select: Adults to be reproduced
% const_struct: a structure to pass constants
% parentnum: the rank of the Adult community

% spike_frac_M: matrix specify the spike fraction 
% rep_num_M: specify the number of Newborns for each spike fraction 
function [comm_rep, H_isolate_out] = pipette_SpikeTest_SPH(comm_selected, newborn_struct, ...
    const_struct, spike_frac_M, rep_num_M, H_isolate_in, parentnum)

comm_rep(1,1:sum(rep_num_M)) = newborn_struct;
gH_max_spike = H_isolate_in.gH_max;
K_HR_spike = H_isolate_in.K_HR;
H_L_spike = H_isolate_in.H_L;
BM_target = const_struct.BM_target;
pcs = const_struct.pcs;

if gH_max_spike<pcs || H_L_spike<pcs
    save('debug_data')
    error ('H isolate is not right')
end

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
rand_idx = randperm(M_counter + H_counter);
all_L_rand = all_L(rand_idx);
all_fp_rand = all_fp(rand_idx);
all_gMax_rand = all_gMax(rand_idx);
all_KR_rand = all_KR(rand_idx);
all_KB_rand = all_KB(rand_idx);

% comm_rep_num = const_struct.comm_rep_num;
% if BM < comm_rep_num*BM_target

if ~isempty(spike_frac_M)
    % calculate the number of cells into each Newborn if the Adult is
    % exhausted
    p_temp = ones(rep_num_M(end), 1) * (BM_target*(1-spike_frac_M(2:end))/BM);
    mn_p = [BM_target*(1-spike_frac_M(1))/BM*ones(1,rep_num_M(1)), transpose(p_temp(:))];
    if sum(mn_p)>1
        save('debug_data')
        error ('mistake in number of Newborns')
    end
    mn_p = [mn_p, 1-sum(mn_p)];
    num_cells = mnrnd(M_counter + H_counter, mn_p);
    num_cells = num_cells(1:end-1);
    H_spike = zeros(sum(rep_num_M), 1);
    temp_counter = 0;
    for s = 1:length(spike_frac_M)
        % assume that the biomass of spiked H cells are 1/log(2)
        H_spike(temp_counter+1:temp_counter+rep_num_M(s)) ...
            = poissrnd(BM_target*spike_frac_M(s)/H_L_spike, [rep_num_M(s), 1]);
        temp_counter = temp_counter + rep_num_M(s);
    end
end
if sum(num_cells) > M_counter+H_counter
    error('not enough cells')
else
    remain_start = sum(num_cells)+1;
    H_idx = find(all_fp_rand(remain_start:end)<0, 1);
    if ~isempty(H_idx)
        H_idx = H_idx+remain_start-1;
        H_isolate_out.H_L = all_L_rand(H_idx);
        H_isolate_out.gH_max = all_gMax_rand(H_idx);
        H_isolate_out.K_HR = all_KR_rand(H_idx);
        if H_isolate_out.gH_max > const_struct.gH_max_Bound+pcs
            save('debug_data')
            error('M isolate is mistaken for H isolate')
        end
    else
%         error('No H isolates')
        H_isolate_out.H_L = -1;
        H_isolate_out.gH_max = -1;
        H_isolate_out.K_HR = -1;
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
    % supplement with H cells with identical biomass and ancestral
    % phenotype
    if H_spike(i)>= 1
        comm_rep(i).H_L(1+H_num : H_spike(i)+H_num) = H_L_spike;
        comm_rep(i).gH_max(1+H_num : H_spike(i)+H_num) = gH_max_spike;
        comm_rep(i).K_HR(1+H_num : H_spike(i)+H_num) = K_HR_spike;
    end
    comm_rep(i).parentnum = parentnum;
end




