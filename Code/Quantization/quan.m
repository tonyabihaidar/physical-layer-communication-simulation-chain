function xq = quan(x, thr, lvl)

%ensure row vector
x = x(:).';
thr = thr(:).';
lvl = lvl(:).';

%ensure correct thr and lvl
if length(lvl) ~= length(thr)+1
    error('Length(lvl) must be length(thr)+1');
end

edges = [-inf, thr, inf];
idx = discretize(x, edges);

xq = lvl(idx);

xq= xq(:).';

end