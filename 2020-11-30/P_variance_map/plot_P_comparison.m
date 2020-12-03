clear
foldername = '3/';
try_num = 3;
ax = [0 1200];
filename = [foldername 'data' num2str(try_num) '.mat'];
load(filename);
P_cal_mat = P_cal*ones(1, 20);
plot(P_cal_mat(:), P_reruns(:), 'ko')
hold on
errorbar(P_cal, mean(P_reruns, 2), std(P_reruns, 1, 2), 'r.', 'markersize', 12, 'linewidth', 2)
plot(ax, ax, 'b:', 'linewidth', 2)
hold off
axis([ax ax])
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'yticklabel',[])
xlabel('P calculated')
ylabel('P simulated')
