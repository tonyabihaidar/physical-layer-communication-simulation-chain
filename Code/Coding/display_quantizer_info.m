function display_quantizer_info(name, xq, thr, lvl, seed)
%DISPLAY_QUANTIZER_INFO  Part 3.1 – Display alphabet, quantizer parameters, and fix seed
% Inputs:
%   name : label for the quantizer (e.g., 'Uniform M=16')
%   xq   : quantized sequence
%   thr  : thresholds
%   lvl  : quantization levels
%   seed : random seed (scalar)

    if nargin < 5, seed = 123; end      
    rng(seed);                     

    M = length(lvl);
    A = lvl;                             % alphabet indices

    fprintf('\n--- %s ---\n', name);
    fprintf('Seed fixed at %d\n', seed);
    fprintf('Alphabet A = {%d}\n', A);
    fprintf('Levels:     '); fprintf('%.4f ', lvl); fprintf('\n');
    fprintf('Thresholds: '); fprintf('%.4f ', thr); fprintf('\n');
    fprintf('Stream length: %d samples\n', numel(xq));
end
