function [xq, thr, lvl, mse] = two_level_quant(x)
%TWO_LEVEL_QUANT  Part (a): Two-level quantization.

    % Use only finite samples
    xf = x(isfinite(x));

    % Median threshold
    thr = median(xf);

    % Compute mean values for the two sides
    lowIdx  = xf <= thr;
    highIdx = xf >  thr;


    % Levels are means of each region
    lvl = [mean(xf(lowIdx)), mean(xf(highIdx))];

    % Use quan function to perform quantization
    xq = quan(x, thr, lvl);

    % Compute MSE
    mse = mean((x(:) - xq(:)).^2, 'omitnan');
end
