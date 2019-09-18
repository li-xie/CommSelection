clear
cycle_num = 1500;
mut_rate = 2e-2;
t_binnum = 340;
mut_tnum = t_binnum*8/10;
mut_event_all = OneCycleReplay_mut(cycle_num, mut_rate/2, 0);
mut_event_res = OneCycleReplay_mut(cycle_num, mut_rate, mut_tnum);
%%
histogram(mut_event_all,(0: max(mut_event_all)));
hold on
histogram(mut_event_res,(0: max(mut_event_res)));
hold off