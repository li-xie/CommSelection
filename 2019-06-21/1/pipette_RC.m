% reproduce Adult communities through pipetting, when six phenotypes could
% mutate

% comm_select: Adults to be reproduced
% const_struct: a structure to pass constants
% dil_factor: dilution factor
% rep_counter: current number of Newborn communities.
% parentnum: the rank of the Adult community
function comm_rep = pipette_RC(comm_selected, comm_struct, const_struct, dil_factor, rep_counter, parentnum)
% maximal number of offspring community from one Adult
comm_rep_num = const_struct.comm_rep_num;
% minimal number of Adults allowed to reproduce
comm_type_num = const_struct.comm_type_num;
% initialize for the Newborn communities
comm_rep(1 : comm_rep_num,1) = comm_struct;
M_L = comm_selected.M_L;
H_L = comm_selected.H_L;
fp = comm_selected.fp;
M_counter = nnz(comm_selected.M_L);
H_counter = nnz(comm_selected.H_L);
c_BM = comm_selected.c_BM;
r_B = comm_selected.r_B;

% assign a random number between 1 and dil_factor to each cell
rand_temp = ceil(rand(M_counter + H_counter,1)*dil_factor);
% if the Adult can generate more Newborns than comm_rep_num, keep only
% comm_rep_numif dil_factor>=comm_rep_num
if dil_factor >= comm_rep_num
    for i = 1:comm_rep_num
        % find all cells assigned with random number i 
        temp_idx = find(rand_temp == i);
        % find M cells
        M_num = nnz(temp_idx <= M_counter);
        % find H cells
        H_num = nnz(temp_idx > M_counter);
        if M_num >= 1
            comm_rep(i).M_L(1:M_num) = M_L(temp_idx(1:M_num));
            comm_rep(i).fp(1:M_num) = fp(temp_idx(1:M_num));
            comm_rep(i).c_BM(1:M_num) = c_BM(temp_idx(1:M_num));
        end
        if H_num >= 1
            comm_rep(i).H_L(1:H_num) = H_L(temp_idx(1+M_num : H_num+M_num)-M_counter);
            comm_rep(i).r_B(1:H_num) = r_B(temp_idx(1+M_num : H_num+M_num)-M_counter);
        end
        comm_rep(i).parentnum = parentnum;
    end
    % keep only a total of comm_type_num*comm_rep_num Newborns
    if rep_counter + comm_rep_num > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    end
else
    for i = 1:dil_factor
        temp_idx = find(rand_temp == i);
        M_num = nnz(temp_idx <= M_counter);
        H_num = nnz(temp_idx > M_counter);
        if M_num >= 1
            comm_rep(i).M_L(1:M_num) = M_L(temp_idx(1:M_num));
            comm_rep(i).fp(1:M_num) = fp(temp_idx(1:M_num));
            comm_rep(i).c_BM(1:M_num) = c_BM(temp_idx(1:M_num));
        end
        if H_num >= 1
            comm_rep(i).H_L(1:H_num) = H_L(temp_idx(1+M_num : H_num+M_num)-M_counter);
            comm_rep(i).r_B(1:H_num) = r_B(temp_idx(1+M_num : H_num+M_num)-M_counter);
        end
        comm_rep(i).parentnum = parentnum;
    end
    if rep_counter + dil_factor > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    else
        comm_rep = comm_rep(1 : dil_factor);
    end
end
    
    
