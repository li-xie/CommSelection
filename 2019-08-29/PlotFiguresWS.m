% 1:3, Neither fixed
% 4:6, N0 fixed
% 7:9, PhiM(0) fixed
% 10:12, both fixed
clear
close all

N=2000; % number of cycles
p0=1;%11.9;

fp=zeros(N,3);
p=zeros(N,3);
% K_MR=zeros(N,3);


TitleName={'Top1Pipette','Top10Pipette','Top1Sort'};
for panel=1:3
    counter=0;
    for i=(panel-1)*3+1:(panel-1)*3+3
        filename=['PlotData/Data' num2str(i) '.mat'];
        load(filename)
        counter=counter+1;
        fp(:,counter)=fp_commmean;
        p(:,counter)=p_commmean/p0;
        %     K_MR(:,counter)=B1_KsFold_commmean;
    end
    
    % plot fp from 3 runs
    figure
    plot((1:N),fp(:,1),'color','k','Linewidth',1);
    hold on
    plot((1:N),fp(:,2),'color','c','Linewidth',1);
    plot((1:N),fp(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
    plot([1 N],[0.41 0.41],'m--','Linewidth',1.5)
    hold off
    axis([1 N 0 0.5])
    set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
    xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
    ylabel('f_P','FontSize',16,'FontName','Arial','fontweight','bold');
    title(TitleName(panel))
    print(['fp-' char(TitleName(panel))],'-dpdf')
    
    % plot P(T) from 3 runs
    figure
    plot((1:N),p(:,1),'color','k','Linewidth',1);
    hold on
    plot((1:N),p(:,2),'color','c','Linewidth',1);
    plot((1:N),p(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
    plot([1 N],[2735.5 2735.5],'m--','Linewidth',1.5)
    hold off
    axis([1 N 0 3000])
    set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
    xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
    ylabel('P(T)','FontSize',16,'FontName','Arial','fontweight','bold');
    title(TitleName(panel))
    print(['PT-' char(TitleName(panel))],'-dpdf')
end
% % plot K_MR from 3 runs
% figure(3)
% plot((1:N),K_MR(:,1),'color','k','Linewidth',1);
% hold on
% plot((1:N),K_MR(:,2),'color','c','Linewidth',1);
% plot((1:N),K_MR(:,3),'color',[0.8 0.8 0.8],'Linewidth',1);
% plot([1 N],[10/3 10/3],'m-','Linewidth',1.5)
% hold off
% axis([1 N 0 5])
% set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
% xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
% ylabel('K_{MR}','FontSize',16,'FontName','Arial','fontweight','bold');