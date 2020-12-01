mult = 1;
% reproduction method
prop = @pipette;
num_cycles = 2;
% number of communities within each cycle
n_tot = 100;
% mutation rate
mu = 2 * 10^-3;
for rep_idx = 1 : 3
    label = sprintf('result_rep%d',rep_idx);
    simulateCommunitySelectionFF(mult,prop,num_cycles,n_tot,mu,label);
end