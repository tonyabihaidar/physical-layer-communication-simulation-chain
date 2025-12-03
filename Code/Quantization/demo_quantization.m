% Demonstration for Part 2 (Quantization)
clear; close all; clc;

% 1. Generate sampled signal
fs = 1000;                 % sampling frequency
t = linspace(0, 1, fs);    % 1 s duration
x = cos(2*pi*5*t);         % clean 5 Hz cosine signal (no noise) 

fprintf('Quantization Demo\n');

% 2. Two-level quantization
[xq2, thr2, lvl2, mse2] = two_level_quant(x);
plot_quant(t, x, xq2, thr2, lvl2, sprintf('Two-level Quantization (MSE=%.6g)', mse2));
fprintf('Two-level quantization:\n');
fprintf('  Threshold: %.6g\n', thr2);
fprintf('  Levels: %.6g  %.6g\n', lvl2);
fprintf('  MSE: %.6g\n\n', mse2);

% 3. Uniform quantization (M = 4, 16, 32)
M_list = [4 16 32];
for M = M_list
    [xqU, thrU, lvlU, mseU] = uniform_quant(x, M);
    plot_quant(t, x, xqU, thrU, lvlU, sprintf('Uniform Quantization M=%d (MSE=%.6g)', M, mseU));
    fprintf('Uniform quantization (M=%d): MSE=%.6g\n', M, mseU);
end

% 4. Lloyd–Max quantization (M = 4, 16, 32)
for M = M_list
    [xqL, thrL, lvlL, mseL, infoL] = lloydmax_quant(x, M);
    plot_quant(t, x, xqL, thrL, lvlL, sprintf('Lloyd–Max Quantization M=%d (MSE=%.6g)', M, mseL));
    fprintf('Lloyd–Max quantization (M=%d): MSE=%.6g, iterations=%d\n', ...
        M, mseL, infoL.iterations);
end

fprintf('\n End of Demo \n');
