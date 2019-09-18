function DataCollect(Fnum,N)
% comm_type_num=10; % number of communities taken after each innoculation cycle
% comm_rep_num=100; % number of replicas of communities
% t_bin=0.05;
% T=17;
% t_binnum=340;
% N=2000;
dataname=['PlotData/Data' num2str(Fnum)];
% beta1c=0.1;
pcs=1e-15;
%%
fp_commmean=zeros(N,1);
% B1_beta1_commstd=zeros(N,1);

% B1_KsFold_commmean=zeros(N,1);
% B1_KsFold_commstd=zeros(N,1);

% B10frac_commmean=zeros(N,1);
% B10frac_commstd=zeros(N,1);
% B1Tfrac_commmean=zeros(N,1);
% B1Tfrac_commstd=zeros(N,1);
% B10_commmean=zeros(N,1);
% B10_commstd=zeros(N,1);
% B1T_commmean=zeros(N,1);
% B1T_commstd=zeros(N,1);
p_commmean=zeros(N,1);
% p_commstd=zeros(N,1);


for n=1:N
% load(['G' num2str(n) '/comm_all.mat']);
load([num2str(Fnum) '/C' num2str(n) '/comm_selected.mat']);
comm_type_num=length(comm_selected);
% B10=zeros(comm_type_num,1);
B1_beta1_mean=zeros(comm_type_num,1);
% B1_beta1_std=zeros(comm_type_num,1);

% B1_KsFold_mean=zeros(comm_type_num,1);
% B1_KsFold_std=zeros(comm_type_num,1);

% B20=zeros(comm_type_num,1);
% B1T=zeros(comm_type_num,1);
% B2T=zeros(comm_type_num,1);

for i=1:comm_type_num
% B10(i)=comm_gen(i).B1_t(1);
% B20(i)=comm_gen(i).B2_t(1);
% B1T(i)=comm_gen(i).B1_t(t_binnum);
% B2T(i)=comm_gen(i).B2_t(t_binnum);
B1_counter=nnz(comm_selected(i).M_L>pcs);
if B1_counter>0
B1_beta1_mean(i)=sum(comm_selected(i).fp.*comm_selected(i).M_L)/sum(comm_selected(i).M_L);
% B1_beta1_std(i)=std(comm_selected(i).B1_beta1(1:B1_counter));

% temp=comm_gen(i).B1_KsFold(1:B1_counter);
% B1_KsFold_mean(i)=sum(comm_gen(i).B1_KsFold.*comm_gen(i).B1_L)/sum(comm_gen(i).B1_L);
% B1_KsFold_std(i)=std(temp((temp<B1_Ks1Max)&(temp>=0.1)));
end

% end
end

fp_commmean(n)=mean(B1_beta1_mean);
% B1_beta1_commstd(n)=std(B1_beta1_mean);
% B1_KsFold_commmean(n)=mean(B1_KsFold_mean);
% B1_KsFold_commstd(n)=std(B1_KsFold_mean);

% B10frac_commmean(n)=mean(B10./(B10+B20));
% B10frac_commstd(n)=std(B10./(B10+B20));
% B1Tfrac_commmean(n)=mean(B1T./(B1T+B2T));
% B1Tfrac_commstd(n)=std(B1T./(B1T+B2T));
% B10_commmean(n)=mean(B10);
% B10_commstd(n)=std(B10);
% B1T_commmean(n)=mean(B1T);
% B1T_commstd(n)=std(B1T);
p_commmean(n)=mean([comm_selected.P]);
% p_commstd(n)=std([comm_selected.p]);
end
% B0_commmean=B10_commmean./B10frac_commmean;

save(dataname,'fp_commmean','p_commmean')