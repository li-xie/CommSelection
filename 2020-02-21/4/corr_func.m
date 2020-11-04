function c = corr_func(x, y)
y_low = mean(y) - 3*std(y);
y_high = mean(y) + 3*std(y);
I = y>y_low & y<y_high;
c = corr(x(I), y(I), 'Type', 'Spearman');