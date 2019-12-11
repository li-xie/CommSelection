clear

spike_frac = 0.3;
comm_type_num = 2;
r = 12;
load('Results.mat');
d = mean(fp0_all(r,:))-mean(fp0_sel(r,:))
load('p30')
load('MfracDiv30')
pp=p';
MfracDivp=MfracDiv';

% fpT_sel=zeros(length(gen2),1);
% B1Tfrac_sel=zeros(length(gen2),1);
% for i=1:length(gen2)
%     fpT_sel(i)=sum(gen2(i).B1_beta1.*gen2(i).B1_L)/sum(gen2(i).B1_L);
%     B1Tfrac_sel(i)=sum(gen2(i).B1_L)/(sum(gen2(i).B1_L)+sum(gen2(i).B2_L));
% end
% range and increment for fp
x0=0.13;
dx=1e-4;
xt=0.16;
% range and increment for B10frac
y0=0.3;
dy=0.005;
yt=0.7;


%%
figure(1)
contour((x0:dx:xt),(y0:dy:yt),pp,(650:50:1.2e3),'linewidth',2)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
% surf((frac_low:1e-2:frac_high),(fp_low:1e-4:fp_high),-Partialmap)
caxis([650 1.2e3])
hold on
% [C, ~] = contour((x0:dx:xt),(y0:dy:yt),MfracDivp,[0 0],'k:','linestyle','none');
load('ContourLine')
plot(C(1, 2:end), C(2, 2:end)*(1-spike_frac),  'color', [1 0.6 0], 'linewidth', 3)
plot(fp0_all(r,:),M0_frac_all(r,:),'ko','linewidth',2)
% plot(fpT_sel,B1Tfrac_sel,'ro','linewidth',2,'markersize',8,'markerfacecolor','r')
plot(fp0_sel(r,:),M0_frac_sel(r,:),'mo','linewidth',2)
plot([mean(fp0_sel(r,:)) mean(fp0_all(r,:))],[0.35 0.35],'o')
% plot([0.1516 0.1516],[y0 yt],'--','color',[0.5 0.5 0.5],'linewidth',2)
hold off
cm = colormap('gray');
c = 1-cm;
colormap(c)
% colorbar
% ylim([min(fp_all) max(fp_all)])
% ylim([0.5 0.9])
axis([x0 xt y0 yt])
xlabel('f_P')
ylabel('phiM')

print(['R' num2str(r) '.pdf'],'-dpdf')
% print('T0=17Colorbar.pdf','-dpdf')
% xlim([0.11 0.135])

