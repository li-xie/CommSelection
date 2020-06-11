clear
c = 2939;
load('C2939/comm_all/newborns')

gH_max_Bound = 0.3;
% the factor for amount of R(0)
V = 1;
% minimal number of Adults allowed to reproduce. comm_type_num = 1 for the
% top-dog strategy, comm_type_num = n for the top n% strategy.
comm_type_num = 10;
% mutation rate corresponding to effective mutation rate of 2e-3
% to turn off the mutation, set mut_rate=0.
mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
Pn_sig = 50; % std of the measurement noise.
N = 100; % number of communities within a cycle
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
T0 = 17;
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4*V; % maximal number of cells in an Adult
nb_popul = BM_target*2; % maximal number of cells in a Newborn
t_bin = 0.05; % time step in the simulation
pcs = 1e-9; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

% ancestral or evolved parameters shown in Table 1
dM = 3.5e-3; % death rate of M
dH = 1.5e-3; % death rate of H
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;

% evolutionary bounds of the phenotypes
gM_max_Bound = 0.7;
fp_Bound = 1; % fp is between 0 and 1
K_MB_Bound = 100/3 ;
K_MR_Bound = 1/3 ;
K_HR_Bound = 0.2 ;


% R(0)=1 for each V
R0 = V;
% if K_MR, K_HR, K_MB is larger than K_singular value, they are effectively
% infinity
K_singular = 1e3;
% structure of an Adult community
comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'gM_max', zeros(max_popul,1), 'K_MB', zeros(max_popul,1), 'K_MR', zeros(max_popul,1),...
    'gH_max', zeros(max_popul,1), 'K_HR', zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of a Newborn community
newborn_struct = struct('M_L',zeros(nb_popul,1),'H_L',zeros(nb_popul,1),'fp',zeros(nb_popul,1),...
    'gM_max', zeros(nb_popul,1), 'K_MB', zeros(nb_popul,1), 'K_MR', zeros(nb_popul,1),...
    'gH_max', zeros(nb_popul,1), 'K_HR', zeros(nb_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of simulation constants
const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target, 't_bin', t_bin,...
    'gM_max_Bound', gM_max_Bound, 'gH_max_Bound', gH_max_Bound,...
    'fp_Bound', fp_Bound, 'K_MB_Bound', K_MB_Bound,'K_MR_Bound', K_MR_Bound, ...
    'K_HR_Bound', K_HR_Bound, 'R0', R0, 'K_singular', K_singular, 'dM', dM, 'dH', dH, ...
    'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH, 'mut_rate', mut_rate);

adults(1:100, 1) = comm_struct;
parfor i = 1:100
    adults(i) = simu_one_comm(newborns(i), comm_struct, const_struct);
end