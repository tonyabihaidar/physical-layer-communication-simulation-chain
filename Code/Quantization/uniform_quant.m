function [xq, thr, lvl, mse, edges] = uniform_quant(x, M)
%UNIFORM_QUANT  Uniform quantizer
%   [xq, thr, lvl, mse, edges] = UNIFORM_QUANT(x, M)
%   - x              : input samples
%   - M              : number of quantization levels
%   - xq             : quantized output (same shape as x)
%   - thr            : thresholds
%   - lvl            : levels
%   - mse            : mean-squared error between x and xq
%   - edges          : edges [xmin ... xmax]

    if nargin < 2, error('Need x and M.'); end
    if M < 2 || M ~= floor(M)
        error('M must be an integer >= 2.');
    end


    xf = x(isfinite(x));
    xmin = min(xf); xmax = max(xf);

    if xmin == xmax
        % Constant signal: all levels collapse to the same value
        edges = repmat(xmin, 1, M+1);
        thr   = edges(2:end-1);          % all equal
        lvl   = repmat(xmin, 1, M);
        xq    = quan(x, thr, lvl);
        mse   = 0;
        return;
    end

    % Uniform finite edges and midpoints (levels)
    edges = linspace(xmin, xmax, M+1);            % M intervals
    lvl   = (edges(1:end-1) + edges(2:end))/2;    % midpoints
    thr   = edges(2:end-1);                        % internal boundaries

    % Quantize using  quan
    xq = quan(x, thr, lvl);

    % MSE
    diff = x(:) - xq(:);
    mse  = mean(diff(~isnan(diff)).^2);

    % Pretty print
    fprintf('M = %d\n', M);
    fprintf('  thresholds (%d): ', numel(thr)); fprintf('%.6g ', thr); fprintf('\n');
    fprintf('  levels     (%d): ', numel(lvl)); fprintf('%.6g ', lvl); fprintf('\n');
    fprintf('  MSE = %.6g\n\n', mse);
end
