% Simulation of community selections 

function simulateCommunitySelectionFF(multiplier,prop_method,num_cycles,...
    num_wells,mu,label)
% Args:
%   multiplier (integer): multiplier for large-volume system; set to 1 for
%       default behavior
%   prop_method (function): propagation method; options are:
%       * @pipette (simulate pipetting; default)
%       * @pipette_const_N (constant newborn population size)
%       * @pipette_const_phiM (constant newborn fraction manufacturer)
%       * @pipette_cNphiM (constant newborn population size and fraction
%         manufacturer)
%   num_cycles (integer): number of cycles to run simulation
%   num_wells (integer): number of wells per cycle
%   mu (float): mutation rate
%   label (string): a label to keep track of this run

if isequal(prop_method,@pipette)
    prop_name = 'pipette';
else
    if isequal(prop_method,@pipette_const_N)
        prop_name = 'const_N';
    else
        if isequal(prop_method,@pipette_const_phiM)
            prop_name = 'const_phiM';
        else
            if isequal(prop_method,@pipette_cNphiM)
                prop_name = 'cNphiM';
            else
                error('you must enter a valid propogation method')
            end
        end
    end
end
str = sprintf('%.0f\n', -log(mu)*10);
resultsfolder = strcat(prop_name,'_m',num2str(multiplier),...
    'w',num2str(num_wells),'pu',str,label,'/');
[~,~,~] = mkdir(resultsfolder(1:end-1));

cycle_duration = 17;
compressWinner = 0; % For now, compression is OFF
digits_fp = 5;
digits_L = 2;
cw_params = [compressWinner, digits_fp, digits_L];
gc.dt = 0.05;
gc.nsteps = ceil(cycle_duration / gc.dt);
gc.b_Hmax = 0.3;
gc.K_HR = 1/5 * multiplier;
gc.K_MR = 1/3 * multiplier;
gc.K_MB = 1/3 * 100 * multiplier;
gc.b_Mmax = 0.7;
gc.c_RH = 10^-4;
gc.c_RM = 10^-4;
gc.c_BM = 1/3;
gc.r_P = 1;
gc.d_M = 3.5 * 10^-3;
gc.d_H = 1.5 * 10^-3;
gc.R_init = 1 * multiplier * 10 ^ 6;
gc.p_mut = mu;
gc.frac_null = 0.5;
gc.sp0 = 0.05;
gc.sn0 = 0.067;
gc.g = 0;
gc.nD = 100;
gc.newborns_per_adult = 10;
if mod(cycle_duration,gc.dt) ~= 0
    error('cycle duration must be a multiple of dt')
end
time = 0 : gc.dt : cycle_duration;
rng('shuffle');
seeds = randi(10^8,[num_cycles,1]);


ic_fp_manu = [0.13;zeros(100 - 1,1)];
ic_L_manu = [1;zeros(100 - 1,1)];
ic_N_manu = [60 * multiplier ;zeros(100 - 1,1)];
ic_L_help = 1;
ic_N_help = 40 * multiplier;
ic_n_genos = 1;
N0 = ic_L_manu' * ic_N_manu + ic_L_help' * ic_N_help;

save_run_conditions(cycle_duration, num_cycles, multiplier, num_wells,...
    cw_params, gc, time, ic_fp_manu, ic_L_manu, ic_N_manu, ic_L_help,...
    ic_N_help, ic_n_genos, resultsfolder, seeds);

newb_fp_manu = cell(1,num_wells);
newb_L_manu = cell(1,num_wells);
newb_N_manu = cell(1,num_wells);
newb_L_help = cell(1,num_wells);
newb_N_help = cell(1,num_wells);
newb_n_genos = cell(1,num_wells);

for i = 1 : num_wells
    newb_fp_manu{i} = ic_fp_manu;
    newb_L_manu{i} = ic_L_manu;
    newb_N_manu{i} = ic_N_manu;
    newb_L_help{i} = ic_L_help;
    newb_N_help{i} = ic_N_help;
    newb_n_genos{i} = ic_n_genos;
end

len_fp_manu_init = length(ic_fp_manu);

wellPlate_fp_manu = {};
wellPlate_L_manu = {};
wellPlate_N_manu = {};
wellPlate_L_help = {};
wellPlate_N_help = {};
wellPlate_Bio_M = {};
wellPlate_Bio_H = {};
wellPlate_R = {};
wellPlate_B = {};
wellPlate_P = {};
wellPlate_n_genos = {};

save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,0,...
    resultsfolder,[]);
for gen = 1 : num_cycles
    % begin maturation step
    rng(seeds(gen));
    tempseeds = randi(10^7,[num_wells,1]);
    parfor i = 1 : num_wells
        rng(tempseeds(i));
        [wellPlate_fp_manu{i},wellPlate_L_manu{i},wellPlate_N_manu{i},...
            wellPlate_L_help{i},wellPlate_N_help{i},wellPlate_Bio_M{i},...
            wellPlate_Bio_H{i},wellPlate_R{i},wellPlate_B{i},wellPlate_P{i},...
            wellPlate_n_genos{i},gc_out]...
            = simulateOneWell(gc,newb_fp_manu{i},newb_L_manu{i},newb_N_manu{i},...
            newb_L_help{i},newb_N_help{i},newb_n_genos{i});
        if ~isequal(gc_out,gc)
            error('You fool! You must never change your global constants')
        end
    end
    % end maturation step; save data
    savedata(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu,...
        wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
        wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos,...
        gen, resultsfolder, cw_params);
    % extract P(T) to identify "winners" for propagation
    P_final = zeros(num_wells,1);
    for i = 1 : num_wells
        P_final(i) = wellPlate_P{i}(end);
    end
    P_sorted = sortrows([P_final,(1:num_wells)'],1);
    % initialize data structures for newborns
    newb_fp_manu = cell(1,num_wells);
    newb_L_manu = cell(1,num_wells);
    newb_N_manu = cell(1,num_wells);
    newb_L_help = cell(1,num_wells);
    newb_N_help = cell(1,num_wells);
    newb_n_genos = cell(1,num_wells);
    % reproduction
    num_wells_filled = 0;
    finished_pipetting = false;
    win_inds = [];
    for i = 1 : num_wells
        windex = round(P_sorted(end - i + 1,2));
        win_inds = [win_inds,windex];
        if cw_params(1)
            [~, lossy_L, lossy_N] = compress(wellPlate_fp_manu{windex},...
                wellPlate_L_manu{windex},wellPlate_N_manu{windex},cw_params(2:3));
            lossy_Bio_M = lossy_L' * lossy_N;
            nD = floor((lossy_Bio_M + wellPlate_Bio_H{windex}(end)) / N0);
        end
        nD = gc.nD; % constant number of dilution recipients
        target_inds = num_wells_filled + 1 : num_wells_filled + gc.newborns_per_adult;
        if target_inds(end) >= num_wells
            target_inds(target_inds > num_wells) = [];
            finished_pipetting = true;
        end
        num_wells_filled = num_wells_filled + gc.newborns_per_adult;
        
        [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos] = prop_method(...
            newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
            wellPlate_fp_manu{windex},wellPlate_L_manu{windex},wellPlate_N_manu{windex},...
            wellPlate_L_help{windex},wellPlate_N_help{windex}, wellPlate_n_genos{windex},...
            target_inds,nD,cw_params,len_fp_manu_init,N0);
        
        if finished_pipetting
            break
        end
    end
    save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,...
        newb_n_genos,gen,resultsfolder,win_inds);
end
end
