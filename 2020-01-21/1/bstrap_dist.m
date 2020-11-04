function dist = bstrap_dist(x,y,f,n_bstraps)
dist = zeros(n_bstraps,1);
for i = 1 : n_bstraps
    idx = randsample(length(x), length(x), true);
    dist(i) = f(x(idx), y(idx));
end
end