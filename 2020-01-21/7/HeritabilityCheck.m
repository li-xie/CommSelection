clear
check_cycle = 2257;
test_rep_max = 50;
q = 0.1;
load(['C' num2str(check_cycle-1) '/ParResults'])
load(['C' num2str(check_cycle) '/OffResults'])
load(['C' num2str(check_cycle-1) '/comm_all/P_all'])
load(['C' num2str(check_cycle-1) '/comm_all/Pn'])
% load(['C' num2str(check_cycle-1) '/comm_all/Pn'])
test_rep_total = min(test_rep_max, size(P_par,2));
f = @(x,y) corr(x, y, 'Type', 'Spearman');
P_par = P_par + Pn_par;
P_off = P_off + Pn_off;

figure(1)
P_sorted = sort(P_all+Pn, 'descend'); 
y = nanmean(P_off(1:3,1:test_rep_total));
null_dist = bstrap_dist(P_sorted(1:test_rep_total)', y', f, 1000);
LB = (quantile(null_dist,q));
UB = (quantile(null_dist,1-q));
ce = corr(P_sorted(1:test_rep_total)', y', 'type', 'spearman');
figure(1)
scatter(P_sorted(1:test_rep_total), mean(P_off(1:3,1:test_rep_total)))
disp([ce, LB, UB])


[P_sorted, I] = sort(P_par(1, :), 'descend');
P_off_mean = nanmean(P_off(4:6,heri_par_idx));
null_dist = bstrap_dist(P_sorted(heri_par_idx)', P_off_mean', f, 1000);
ce = f(P_sorted(heri_par_idx)', P_off_mean');
LB = (quantile(null_dist,q));
UB = (quantile(null_dist,1-q));
disp([ce, LB, UB])
figure(2)
scatter(P_sorted(1:test_rep_total), P_off_mean)


[P_sorted, I] = sort(P_par(2, :), 'descend');
P_off_mean = nanmean(P_off(7:9,:));
null_dist = bstrap_dist(P_sorted(1:test_rep_total)', P_off_mean', f, 1000);
ce = f(P_sorted(1:test_rep_total)', P_off_mean');
LB = (quantile(null_dist,q));
UB = (quantile(null_dist,1-q));
disp([ce, LB, UB])
figure(3)
scatter(P_sorted(1:test_rep_total), P_off_mean)

function dist = bstrap_dist(x,y,f,n_bstraps)
dist = zeros(n_bstraps,1);
for i = 1 : n_bstraps
    idx = randsample(length(x), length(x), true);
    dist(i) = f(x(idx), y(idx));
end
end