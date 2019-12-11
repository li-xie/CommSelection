clear
comm_type_num = 2;
comm_rep_num = 50;
rp = 20;
fp0_all = zeros(rp, comm_type_num*comm_rep_num);
fp0_sel = zeros(rp, comm_type_num);
P_all = zeros(rp, comm_type_num*comm_rep_num);
M0_frac_all = zeros(rp, comm_type_num*comm_rep_num);
M0_frac_sel = zeros(rp, comm_type_num);
for j = 1:rp
    folder_name = ['R' num2str(j)];
    load([folder_name '/newborns'])
    load([folder_name '/adults'])
    for i=1:comm_type_num*comm_rep_num
        fp0_all(j,i) = sum(newborns(i).fp.*newborns(i).M_L)/sum(newborns(i).M_L);
        M0_frac_all(j,i) = sum(newborns(i).M_L)/(sum(newborns(i).M_L)+sum(newborns(i).H_L));
    end
    P_all(j, :) = [adults.P];
    [~, idx] = sort(P_all(j, :), 'descend');
    fp0_sel(j, :) = fp0_all(j, idx(1:comm_type_num));
    M0_frac_sel(j, :) = M0_frac_all(j, idx(1:comm_type_num));
end
save('Results','fp0_all','fp0_sel','M0_frac_all','M0_frac_sel','P_all')    
    