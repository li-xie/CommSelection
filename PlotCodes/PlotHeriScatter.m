clear
if ~exist('SVG/1Figs','dir')
    mkdir('SVG/1Figs')
end
p_flag = true;
fnum = 1;
check_cycle = 630;
check_time = 4;
test_rep_max = 100;
off_rep_max = 6;
ax = [50 550 50 550];
q = 0.05;
load(['PlotData/Check' num2str(fnum) '/C' num2str(check_cycle) '/ParResults'])
load(['PlotData/Check' num2str(fnum) '/C' num2str(check_cycle) '/OffResults'])
load(['PlotData/Check' num2str(fnum) '/C' num2str(check_cycle) '/P_all'])
load(['PlotData/Check' num2str(fnum) '/C' num2str(check_cycle) '/Pn'])
load(['PlotData/Check' num2str(fnum) '/CheckSummary.mat'])
test_rep_total = min(test_rep_max, size(P_par,2));
f = @(x,y) corr(x, y, 'Type', 'Spearman');
P_par = P_par + Pn_par;
P_off = P_off + Pn_off;

figure(1)
P_sorted = sort(P_all+Pn, 'descend'); 
P_off_mean = nanmedian(P_off(1:off_rep_max,1:test_rep_total));
% ce = f(P_sorted(heri_par_idx)', P_off_mean')
figure(1)
scatter(P_sorted(1:test_rep_total), nanmedian(P_off(1:off_rep_max,1:test_rep_total)), 'k', 'filled')
axis(ax);
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
% disp([corr_mean(check_time, 1), (ub_m(check_time, 1)-lb_m(check_time, 1))/2])
if p_flag == 1
    print('SVG/1Figs/Heri0.svg', '-dsvg')
end
    
for i = 1:2
    [P_sorted, I] = sort(P_par(i, :), 'descend');
    P_off_mean = nanmedian(P_off(off_rep_max*i+1 : off_rep_max*(i+1), heri_par_idx));
%     ce = f(P_sorted(heri_par_idx)', P_off_mean')
    figure(i+1)
    scatter(P_sorted(heri_par_idx), P_off_mean, 'k', 'filled')
    axis(ax);
    set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
    if p_flag == 1
        print(['SVG/1Figs/Heri' num2str(i) '.svg'], '-dsvg')
    end
end

function dist = bstrap_dist(x,y,f,n_bstraps)
dist = zeros(n_bstraps,1);
for i = 1 : n_bstraps
    idx = randsample(length(x), length(x), true);
    dist(i) = f(x(idx), y(idx));
end
end