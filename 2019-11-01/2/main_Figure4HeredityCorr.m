% This code generates Fig 4 B_D by reproducing every Adult and examining
% the average fp and phi_M(0) of their offsprings
clear
% reproduce every Adults from the cycle of cycle_num
cycle_num = 180;
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
N = 100;
T0 = 17; % maturation time
comm_type_num = 1; % minimal number of Adults allowed to reproduce
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4; % maximal number of cells in the community
t_bin = 0.05; % time step in the simulation
pcs=1e-15; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target);
% initialize matrice for the following data:
% average fp(0) of parents 
fp0_p = zeros(comm_type_num*comm_rep_num, 1); 
% average fp(0) of offspring
fp0_o = zeros(comm_type_num*comm_rep_num, comm_type_num*comm_rep_num);
% phi_M(0) of parents 
phi_M0_p = zeros(comm_type_num*comm_rep_num, 1);
% phi_M(T) of parents 
phi_MT_p = zeros(comm_type_num*comm_rep_num, 1);
% phi_M(0) of offspring
phi_M0_o = zeros(comm_type_num*comm_rep_num, comm_type_num*comm_rep_num);
% BM(0) of parents
BM0_p = zeros(comm_type_num*comm_rep_num, 1);
% BM(0) of offsprings
BM0_o = zeros(comm_type_num*comm_rep_num, comm_type_num*comm_rep_num);
% mean and std of fp(0) among offspring from each Adult
fp0_os = zeros(2, comm_type_num*comm_rep_num);
% mean and std of phi_M(0) among offspring from each Adult
phi_M0_os=zeros(2, comm_type_num*comm_rep_num);
% mean and std of BM(0) among offspring from each Adult
BM0_os = zeros(2, comm_type_num*comm_rep_num);

foldername = ['C' num2str(cycle_num)];
load([foldername '/distrng.mat']);
foldername = ['C' num2str(cycle_num) '/comm_all'];
load([foldername '/adults.mat']);
load([foldername '/newborns.mat']);
% set the random number generator state
rng(distrng);

% reproduce every Adults 
[~, I] = sort([adults.P], 'descend');
for j = 1:comm_type_num*comm_rep_num
    comm_temp=adults( I(j) );
    dil_factor = floor( (sum(comm_temp.M_L)+sum(comm_temp.H_L))/BM_target );
    if dil_factor==0
        continue
    end
    % limit the number of Newborns to the smaller of dil_factor and comm_rep_num
    rep_num_temp = min(dil_factor, comm_rep_num);
    newborns_temp = pipette(comm_temp, comm_struct, const_struct, dil_factor, 0, j);
    % calculate the average fp(0), BM(0) and phi_M(0) of all the Newborns
    % from an Adult
    fp0_temp = zeros(comm_type_num*comm_rep_num, 1);
    BM0_temp = zeros(comm_type_num*comm_rep_num, 1);
    phi_M0_temp = zeros(comm_type_num*comm_rep_num, 1);
    for k = 1 : rep_num_temp
        fp0_temp(k) = sum(newborns_temp(k).fp.*newborns_temp(k).M_L)/sum(newborns_temp(k).M_L);
        BM0_temp(k) = sum(newborns_temp(k).M_L+newborns_temp(k).H_L);
        phi_M0_temp(k) = sum(newborns_temp(k).M_L)/BM0_temp(k);
    end
    fp0_o(:,j) = fp0_temp;
    phi_M0_o(:,j) = phi_M0_temp;
    BM0_o(:,j) = BM0_temp;
end
%%
% number of offspring from each Adult
offspring_num = zeros(comm_type_num*comm_rep_num, 1);
% calculate the average fp(0), BM(0) and phi_M(0) of all the parent
% communities
for j=1:comm_type_num*comm_rep_num
    fp0_p(j) = sum(newborns(I(j)).fp.*newborns(I(j)).M_L)/sum(newborns(I(j)).M_L);
    BM0_p(j) = sum(newborns(I(j)).M_L+newborns(I(j)).H_L);
    phi_M0_p(j) = sum(newborns(I(j)).M_L)/BM0_p(j);
    phi_MT_p(j) = sum(adults(I(j)).M_L)/(sum(adults(I(j)).M_L)+sum(adults(I(j)).H_L));
    offspring_idx = BM0_o(:, j)>pcs;
    offspring_num(j) = nnz(offspring_idx);
    fp0_os(1, j) = mean(fp0_o(offspring_idx,j));
    phi_M0_os(1, j) = mean(phi_M0_o(offspring_idx,j));
    BM0_os(1, j) = mean(BM0_o(offspring_idx,j));
    fp0_os(2, j) = std(fp0_o(offspring_idx,j));
    phi_M0_os(2, j) = std(phi_M0_o(offspring_idx,j));
    BM0_os(2, j) = std(BM0_o(offspring_idx,j));
end
%%
close all
figure(1)
plot(fp0_p, fp0_os(1,:), 'ko', 'linewidth', 2, 'markersize', 6)
axis([0.14 0.15 0.13 0.15])
set(gca, 'LineWidth', 2, 'FontSize', 16, 'FontName', 'Arial',...
    'fontweight', 'bold', 'units', 'inches', ...
    'position', [1 1 3 3], 'ticklength', [0.04 0.04])
xlabel('fp(0) of previous cycle')
ylabel('mean fp(0) of current cycle')
title('Fig4B')

figure(2)
plot(BM0_p,BM0_os(1,:),'ko','linewidth',2,'markersize',6)
axis([50 150 75 125])
set(gca, 'LineWidth', 2, 'FontSize', 16, 'FontName', 'Arial',...
    'fontweight', 'bold', 'units', 'inches', ...
    'position', [1 1 3 3], 'ticklength', [0.04 0.04])
xlabel('BM(0) of previous cycle')
ylabel('mean BM(0) of current cycle')
title('Fig4C')

figure(3)
plot(phi_M0_p, phi_M0_os(1,:), 'ko', 'linewidth', 2, 'markersize', 6)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
axis([0.5 0.9 0.6 0.8])
set(gca, 'LineWidth', 2, 'FontSize', 16, 'FontName', 'Arial',...
    'fontweight', 'bold', 'units', 'inches', ...
    'position', [1 1 3 3], 'ticklength', [0.04 0.04])
xlabel('phi_M(0) of previous cycle')
ylabel('mean phi_M(0) of current cycle')
title('Fig4D')
        
