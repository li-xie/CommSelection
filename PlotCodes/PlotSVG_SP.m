SimuCond = {'5HM_Heri', '5HM_CT', '5HM_Rand'};
cut_file_path = 'svg_cut.svg';
f_num = [1:3; 4:6; 7:9];
N = 3000;
print_flag = true;
spike_flag = true;
for i = 1:length(SimuCond)
    PlotSP_func(N, f_num(i,:), SimuCond{i}, print_flag, spike_flag)
    svg_name = dir(SimuCond{i});
    svg_name = svg_name(3:end);
    s_num = length(svg_name);
    for j = 1:s_num
        svg_bgRM([SimuCond{i} '/' svg_name(j).name], cut_file_path)
    end
end