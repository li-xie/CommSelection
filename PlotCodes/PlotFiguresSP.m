clear
% close all

N = 3000; % number of cycles
N1 = 200; % number of cycles to zoom in
print_flag = true; % whether or not print the figures.
rc = true;
f_num = [1:3];
pref = 'Heri5HM';

if ~exist(pref, 'dir')
    mkdir(pref)
end
fp = zeros(N,3);
P = zeros(N,3);
K_MR = zeros(N,3);
K_MB = zeros(N,3);
K_HR = zeros(N,3);
b_Mmax = zeros(N,3);
b_Hmax = zeros(N,3);
B1Tfrac = zeros(N,3);

% TitleName={'Neither fixed','N0 fixed','PhiM(0) fixed','Both fixed'};
counter = 0;
for i = 1:3
    %(panel-1)*3+1:(panel-1)*3+3
    filename=['PlotData/Data' num2str(f_num(i)) '.mat'];
    %     filename=['PlotData/Data' num2str(i) '.mat'];
    load(filename)
    counter=counter+1;
    fp(:,counter) = fp_commmean;
    P(:,counter) = p_commmean;
    K_MR(:,counter) = K_MR_commmean;
    K_MB(:,counter) = K_MB_commmean;
    K_HR(:,counter) = K_HR_commmean;
    b_Mmax(:,counter) = g_Mmax_commmean;
    b_Hmax(:,counter)=g_Hmax_commmean;
end
fc = 0; % counter for figures
%%
%     fc = fc+1;
%     figure(fc)
% plot((1:N),B1Tfrac(:,1),'color','k','Linewidth',1);
% hold on
% plot((1:N),B1Tfrac(:,2),'color','c','Linewidth',1);
% plot((1:N),B1Tfrac(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
% % plot([1 N],[0.5 0.5],'g--','Linewidth',1.5)
% axis([0 1500 0 1])
% hold off
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
% xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('M fraction','FontSize',16,'FontName','Arial','fontweight','bold');
% % title(TitleName(panel))
% print('B1frac(T)','-dpdf')
%%
% plot P(T) from 3 runs
fc = fc+1;
figure(fc)
plot((1:N),P(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),P(:,2),'color','c','Linewidth',1);
plot((1:N),P(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[2735.5 2735.5],'m--','Linewidth',1.5)
hold off
axis([1 N 0 3000])
xticks((1e3:1e3:N))
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('P(T)','FontSize',16,'FontName','Arial','fontweight','bold');
% title(TitleName(panel))
if print_flag == 1
    print([pref '/' 'P(T).svg'],'-dsvg')
end
%%
% plot fp from 3 runs
fc = fc+1;
figure(fc)
plot((1:N),fp(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),fp(:,2),'color','c','Linewidth',1);
plot((1:N),fp(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[0.41 0.41],'m--','Linewidth',1.5)
hold off
xticks((1e3:1e3:N))
axis([1 N 0 0.5])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('f_P','FontSize',16,'FontName','Arial','fontweight','bold');
% title(TitleName(panel))
if print_flag == 1
    print([pref '/' 'fp.svg'],'-dsvg')
end
%%
fc = fc+1;
figure(fc)
plot((1:N),b_Mmax(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),b_Mmax(:,2),'color','c','Linewidth',1);
plot((1:N),b_Mmax(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot((1:N),b_Hmax(:,1),'color','k','Linewidth',1);
plot((1:N),b_Hmax(:,2),'color','c','Linewidth',1);
plot((1:N),b_Hmax(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[0.7 0.7],'g--','Linewidth',1.5)
plot([1 N],[0.3 0.3],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:N))
axis([1 N 0 0.85])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('b_{M,Hmax}','FontSize',16,'FontName','Arial','fontweight','bold');
if print_flag == 1
    print([pref '/' 'gMax_MH.svg'],'-dsvg')
end
% % plot b_Mmax from 3 runs
%  fc = fc+1;
%  figure(fc)
% plot((1:N1),b_Mmax(1:N1,1),'color','k','Linewidth',1);
% hold on
% plot((1:N1),b_Mmax(1:N1,2),'color','c','Linewidth',1);
% plot((1:N1),b_Mmax(1:N1,3),'color',[0.8 0.8 0.8],'Linewidth',1);
% plot([1 N1],[0.7 0.7],'g--','Linewidth',1.5)
% hold off
% axis([1 N1 0 0.85])
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
% xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('b_{Mmax}','FontSize',16,'FontName','Arial','fontweight','bold');
% % title(TitleName(panel))
% if print_flag == 1
%    print([pref '/' 'gMax_M.svg'],'-dsvg')
% end

% %%
% % plot b_Hmax from 3 runs
%  fc = fc+1;     figure(fc)
% plot((1:N1),b_Hmax(1:N1,1),'color','k','Linewidth',1);
% hold on
% plot((1:N1),b_Hmax(1:N1,2),'color','c','Linewidth',1);
% plot((1:N1),b_Hmax(1:N1,3),'color',[0.8 0.8 0.8],'Linewidth',1);
% plot([1 N1],[0.8 0.8],'g--','Linewidth',1.5)
% hold off
% axis([1 N1 0 0.85])
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
% xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('b_{Hmax}','FontSize',16,'FontName','Arial','fontweight','bold');
% % title(TitleName(panel))
% if print_flag == 1
%    print([pref '/' 'gMax_M.svg'],'-dsvg')
% end
%%
% plot K_MR from 3 runs
fc = fc+1;
figure(fc)
plot((1:N),1./K_MR(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),1./K_MR(:,2),'color','c','Linewidth',1);
plot((1:N),1./K_MR(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[3 3],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:N))
axis([1 N 0 4])
% yticks([0 0.1 0.2 0.3 0.4])
% yticklabels({'0','1','2','3','4' })
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('1/K_{MR}','FontSize',16,'FontName','Arial','fontweight','bold');
% title(TitleName(panel))
if print_flag == 1
    print([pref '/' 'K_MR-1.svg'],'-dsvg')
end
% plot K_MB from 3 runs
fc = fc+1;
figure(fc)
plot((1:N),1./K_MB(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),1./K_MB(:,2),'color','c','Linewidth',1);
plot((1:N),1./K_MB(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[0.030 0.030],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:N))
yticks([0 0.006 0.012 0.018 0.024 0.03])
yticklabels({'0','6','12','18','24','30'})
axis([1 N 0 0.033])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('1/K_{MB}','FontSize',16,'FontName','Arial','fontweight','bold');
% title(TitleName(panel))
if print_flag == 1
    print([pref '/' 'K_MB-1.svg'],'-dsvg')
end

% plot K_HR from 3 runs
fc = fc+1;
figure(fc)
plot((1:N),1./K_HR(:,1),'color','k','Linewidth',1);
hold on
plot((1:N),1./K_HR(:,2),'color','c','Linewidth',1);
plot((1:N),1./K_HR(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
plot([1 N],[5 5],'g--','Linewidth',1.5)
hold off
axis([1 N 0 6])
xticks((1e3:1e3:N))
% yticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
% yticklabels({'0','1','2','3','4','5','6' })
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('1/K_{HR}','FontSize',16,'FontName','Arial','fontweight','bold');
% title(TitleName(panel))
if print_flag == 1
    print([pref '/' 'K_HR-1.svg'],'-dsvg')
end
%%
if rc == true
    c = {'k', 'c', [0.8 0.8 0.8]};
    fc = fc+1;
    figure(fc)
    for j = 1:3
        load(['PlotData/Check' num2str(f_num(j)) '/CheckSummary'])
        frac = zeros(N, 1);
        check_cycle_m = [0 check_cycle_m];
        for i = 1:length(check_cycle_m)-1
            frac(check_cycle_m(i)+1:check_cycle_m(i+1)) = spike_before_m(i, 1);
        end
        frac(check_cycle_m(end):N) = spike_after_m(end, 1);
        plot((1:N), frac, 'color', c{j})
        hold on
    end
    hold off
    xticks((1e3:1e3:N))
    yticks([-0.6 -0.3 0 0.3 0.6])
    axis([1 N -0.7 0.7])
    set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04],'yticklabel',[])
    xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
    %     ylabel('Spikeing fraction','FontSize',16,'FontName','Arial','fontweight','bold');
    if print_flag == 1
        print([pref '/' 'SpikeRatio.svg'],'-dsvg')
    end
end