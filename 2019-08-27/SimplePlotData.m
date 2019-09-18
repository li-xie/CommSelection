clear
% close all
figure(1)
for i=7:9
    load(['PlotData/Data' num2str(i)])
    plot(fp_commmean,'.-')
    hold on
end
hold off

% figure(2)
% for i=10:12
%     load(['PlotData/Data' num2str(i)])
%     plot(B1_beta1_commmean,'.-')
%     hold on
% end
% hold off