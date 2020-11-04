% trim zeros from the Newborns structure
function newborns = newborn_trim(comm_all,BM_target,pcs)
comm_num = length(comm_all);
newborns(1 : comm_num,1)=struct('M_L',zeros(BM_target, 1), 'H_L', zeros(BM_target,1),...
    'fp', zeros(BM_target, 1), 'c_BM', zeros(BM_target,1),'r_B', zeros(BM_target,1),...
    'parentnum', 0, 'rseed', uint32(0));
for i = 1:comm_num
    cell_num = nnz(comm_all(i).M_L > pcs);
    newborns(i).M_L(1:cell_num) = comm_all(i).M_L(1:cell_num);
    newborns(i).fp(1:cell_num) = comm_all(i).fp(1:cell_num);
    newborns(i).c_BM(1:cell_num) = comm_all(i).c_BM(1:cell_num);
    cell_num = nnz(comm_all(i).H_L > pcs);
    newborns(i).H_L(1:cell_num) = comm_all(i).H_L(1:cell_num);
    newborns(i).r_B(1:cell_num) = comm_all(i).r_B(1:cell_num);
        newborns(i).parentnum = comm_all(i).parentnum;
    newborns(i).rseed = comm_all(i).rseed;
end