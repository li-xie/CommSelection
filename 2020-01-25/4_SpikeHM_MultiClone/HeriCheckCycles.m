clear
load('comm_all/check_cycle_m')
check_counter = length(check_cycle_m);
spike_before_m = zeros(check_counter, 5);
spike_after_m = zeros(check_counter, 5);
lb_m = zeros(check_counter, 5);
ub_m = zeros(check_counter, 5);
corr_mean = zeros(check_counter, 5);
for i = 1:check_counter
    n = check_cycle_m(i);
    load(['C' num2str(n-1) '/ParResults'])
    load(['C' num2str(n) '/OffResults'])
    spike_before_m(i, :) = spike_all;
    corr_mean(i, :) = heri';
    lb_m(i, :) = lb';
    ub_m(i, :) = ub';
    load(['C' num2str(n) '/spike_all'])
    spike_after_m(i, :) = spike_all;
end
%%
plot(check_cycle_m, spike_before_m(:,1), '.-')
    
    