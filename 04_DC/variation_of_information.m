function sigma = variation_of_information(X,Y)
uniX = unique(X)';
uniY = unique(Y)';
n = numel(X);
sigma = 0;
for ix = uniX
    x = X == ix;
    p = sum(x)/n;
    for iy = uniY
        y = Y == iy;
        q = sum(y)/n;
        r = sum(x & y)/n;
        if r > 0
            sigma = sigma + (r * (log2(r/p) + log2(r/q)));
        end
    end
end
sigma = abs(sigma);
end