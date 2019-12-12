% trim zeros from the Newborns structure
function newborns = newborn_trim(comm_all,BM_target,pcs)
comm_num = length(comm_all);
newborns(1 : comm_num,1)=struct('M_L',zeros(BM_target, 1), 'H_L', zeros(BM_target,1),...
    'fp', zeros(BM_target, 1), 'gM_max',zeros(BM_target, 1),'K_MB', zeros(BM_target,1),...
    'K_MR', zeros(BM_target, 1), 'gH_max', zeros(BM_target, 1),'K_HR', zeros(BM_target,1),...
    'parentnum', 0, 'rseed', uint32(0));
for i = 1:comm_num
    cell_num = nnz(comm_all(i).M_L > pcs);
    newborns(i).M_L(1:cell_num) = comm_all(i).M_L(1:cell_num);
    newborns(i).fp(1:cell_num) = comm_all(i).fp(1:cell_num);
    newborns(i).gM_max(1:cell_num) = comm_all(i).gM_max(1:cell_num);
    newborns(i).K_MB(1:cell_num) = comm_all(i).K_MB(1:cell_num);
    newborns(i).K_MR(1:cell_num) = comm_all(i).K_MR(1:cell_num);
    cell_num = nnz(comm_all(i).H_L > pcs);
    newborns(i).H_L(1:cell_num) = comm_all(i).H_L(1:cell_num);
    newborns(i).gH_max(1:cell_num) = comm_all(i).gH_max(1:cell_num);
    newborns(i).K_HR(1:cell_num) = comm_all(i).K_HR(1:cell_num);
    newborns(i).parentnum = comm_all(i).parentnum;
    newborns(i).rseed = comm_all(i).rseed;
end