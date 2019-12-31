% clear
check_cycle = 1700;
load(['C' num2str(check_cycle-1) '/ParResults'])
load(['C' num2str(check_cycle) '/OffResults'])
load(['C' num2str(check_cycle-1) '/comm_all/P_all'])
load(['C' num2str(check_cycle-1) '/comm_all/Pn'])
figure(1)
[~, I] = sort(P_all+Pn, 'descend');
P_sorted = P_all(I);
figure(1)
scatter(P_sorted, mean(P_off(1:3,:)))
temp = corrcoef(P_sorted, mean(P_off(1:3,:)));
title(num2str(temp(1,2)))
figure(2)
scatter(P_par(1,:), mean(P_off(4:6,:)))
temp = corrcoef(P_par(1,:), mean(P_off(4:6,:)));
title(num2str(temp(1,2)))
figure(3)
scatter(P_par(2,:), mean(P_off(7:9,:)))
temp = corrcoef(P_par(2,:), mean(P_off(7:9,:)));
title(num2str(temp(1,2)))