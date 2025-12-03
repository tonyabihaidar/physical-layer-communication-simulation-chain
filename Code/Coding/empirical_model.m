function [p,H] = empirical_model(name, xq, lvl)
% (b) Empirical modeling: p(a) and H(A)

    [tf, idx] = ismember(xq, lvl);
    if ~all(tf), error('Some samples not in lvl.'); end

    M = numel(lvl);
    N = numel(idx);
    counts = zeros(M,1);
    for k = 1:length(idx)
        counts(idx(k)) = counts(idx(k)) + 1;
    end
    p = counts / N;

    H = -sum(p(p>0) .* log2(p(p>0)));

    fprintf('\n--- %s: Empirical model ---\n', name);
    fprintf('Level\tp(level)\n');
    for i = 1:M
        fprintf('%8.4f\t%.4f\n', lvl(i), p(i));
    end
    fprintf('Entropy H(A) = %.4f bits/symbol\n', H);
end
