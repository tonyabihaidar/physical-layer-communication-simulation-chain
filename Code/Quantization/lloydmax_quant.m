function [xq, thr, lvl, mse] = lloydmax_quant(x, M, maxIter, tol)
%LLOYDMAX_QUANT  Lloyd–Max quantizer learned from empirical samples.
%   [xq, thr, lvl, mse, info] = LLOYDMAX_QUANT(x, M, maxIter, tol)
%   Inputs:
%     x       : input sample
%     M       : number of quantization levels
%     maxIter : max iterations, default 100
%     tol     : convergence tolerance for levels/thresholds, default 1e-9
%
%   Outputs:
%     xq   : quantized signal
%     thr  : thresholds
%     lvl  : levels
%     mse  : mean-squared error
%     info : struct with fields:
%              .iterations  : number of iterations executed
%              .converged   : logical
%              .emptyBins   : indices of bins that were empty at least once
%              .init        : 'uniform'

    if nargin < 2, error('Need x and M.'); end
    if M < 2 || M ~= floor(M), error('M must be an integer >= 2.'); end
    if nargin < 3 || isempty(maxIter), maxIter = 100; end
    if nargin < 4 || isempty(tol),     tol     = 1e-9; end

    % Flatten to a vector for training; keep original shape to restore later
    xv = x(:);
    mask = isfinite(xv);
    xtrain = xv(mask);

    xmin = min(xtrain); xmax = max(xtrain);
    if xmin == xmax
        % Constant signal → trivial quantizer
        lvl = repmat(xmin, 1, M);
        thr = repmat(xmin, 1, M-1);   % any set of equal thresholds
        xq  = quan(x, thr, lvl);
        mse = 0;
        return;
    end

    % Initialization (uniform over data range)
    edges = linspace(xmin, xmax, M+1);
    lvl   = (edges(1:end-1) + edges(2:end))/2;   % midpoints
    thr   = (lvl(1:end-1) + lvl(2:end))/2;       % initial thresholds

    % Lloyd–Max iterations
    for it = 1:maxIter
        % Assign samples to regions using current thresholds
        % Regions are (-Inf, thr(1)], (thr(1), thr(2)], ..., (thr(M-1), Inf)
        idx = discretize(xtrain, [-inf, thr, inf]);   % 1..M

        % Update levels as conditional means
        lvl_new = lvl;  % copy for potential empty-bin handling
        emptyNow = [];

        for k = 1:M
            sel = (idx == k);
            if any(sel)
                lvl_new(k) = mean(xtrain(sel));
            else
                % Empty bin: fix by splitting the largest-populated neighbor's span midpoint
                emptyNow(end+1) = k; %#ok<AGROW>
                % Simple fallback: keep previous level (stable) to avoid NaN,
                % or, if possible, interpolate between neighbors:
                if k > 1 && k < M
                    lvl_new(k) = 0.5*(lvl(k-1) + lvl(k+1));
                elseif k == 1
                    lvl_new(k) = lvl(k) - 0.1*(xmax - xmin);
                else % k == M
                    lvl_new(k) = lvl(k) + 0.1*(xmax - xmin);
                end
            end
        end

        % Update thresholds as midpoints between adjacent levels
        thr_new = 0.5*(lvl_new(1:end-1) + lvl_new(2:end));

        % Check convergence
        dLvl = max(abs(lvl_new - lvl));
        dThr = max(abs(thr_new - thr));
        lvl = lvl_new;
        thr = thr_new;

        if max(dLvl, dThr) < tol
            break;
        end
    end

    % Quantize with the learned quantizer
    xq = quan(x, thr, lvl);

    % MSE 
    diff = x(:) - xq(:);
    mse  = mean(diff(~isnan(diff)).^2);
end
