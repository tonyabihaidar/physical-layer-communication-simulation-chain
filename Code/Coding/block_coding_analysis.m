function [L_fixed_symbol,avg_len_per_symbol_block,Hk_rate]=block_coding_analysis(name, xq, lvl, k, H_symbol, L_huff_symbol)
% Inputs:
%   name          : Label for the quantizer (e.g., 'Uniform M=16')
%   xq            : Quantized sequence (as levels)
%   lvl           : Quantizer levels
%   k             : Block size (e.g., 2, 3)
%   H_symbol      : Single-symbol entropy H(A) from empirical_model.m
%   L_huff_symbol : Avg. Huffman length (symbol-wise) from huffman_encode.m

M = numel(lvl);
N = numel(xq);

% Map original sequence of levels to a sequence of indices (1..M)
[tf, idx] = ismember(xq, lvl);
if ~all(tf)
    error('block_coding_analysis: some samples in xq are not present in lvl.');
end

% ---- Guard: k must be a positive integer
if ~(isscalar(k) && k == floor(k) && k > 0)
    error('k must be a positive integer.');
end

%part e)
num_blocks = floor(N / k);
idx_truncated = idx(1 : num_blocks * k); %this is used to discard any symbols
%that do not form a full block

% ---- Guard: no full blocks case (N < k)
if num_blocks == 0
    warning('block_coding_analysis: N < k (no full blocks). Returning NaNs.');
    L_fixed_symbol = ceil(log2(M));
    avg_len_per_symbol_block = NaN;
    Hk_rate = NaN;
    return;
end

block_matrix = reshape(idx_truncated, k, num_blocks);
%now we have a k x num_block matrix where each column is a block

%now to find the unique blocks, there's a built in function to do so in
%matrix form but it only works row wise so we'll have to transpose using
%the ' command which transposes any given matrix    
[unique_blocks_T, ~, block_indices] = unique(block_matrix', 'rows');
num_unique_blocks = size(unique_blocks_T, 1);

%now to calculate the probability of each unique code
counts = histcounts(block_indices, 1:num_unique_blocks+1);
p_blocks = counts' / num_blocks;

%part f%)

L_fixed_symbol = ceil(log2(M));
fprintf('Fixed-Length (Symbol-wise)  \t%.4f\n', L_fixed_symbol);

% ---- Added (kept separate from returns): Fixed-length baseline for BLOCKS
% bits per block = ceil(log2 K), per-symbol baseline = that / k
L_fixed_block          = ceil(log2(max(1,num_unique_blocks)));   % bits per block
L_fixed_block_per_sym  = L_fixed_block / k;                      % per symbol
fprintf('Fixed-Length (Block-wise)/k \t%.4f\n', L_fixed_block_per_sym);

% 5. Huffman Code (Block-wise)
%    Design a new Huffman code for the alphabet of blocks.
%    The "symbols" are integers from 1 to num_unique_blocks.
dict_block = huffmandict(num2cell(1:num_unique_blocks), p_blocks);
%the following built in MATLAB funciton calculates how many bits in each
%codeword
codeword_lengths = cellfun(@length, {dict_block{:,2}});
% Calculate average length per block
avg_len_block = sum(p_blocks .* codeword_lengths');
% Normalize to get average length per original symbol
avg_len_per_symbol_block = avg_len_block / k; %Lk

% Calculate the block entropy (Hk) and the entropy rate (Hk/k)
mask = p_blocks > 0;
Hk = -sum(p_blocks(mask) .* log2(p_blocks(mask)));
Hk_rate = Hk / k;

fprintf('\n--- %s: Summary for part (g) ---\n', name);
fprintf('Entropy H(A) (per symbol)                          : %.4f bits/symbol\n', H_symbol);
fprintf('Average length using Huffman (symbol-wise)         : %.4f bits/symbol\n', L_huff_symbol);
fprintf('Block entropy rate H_k/k (per symbol, k=%d)        : %.4f bits/symbol\n', k, Hk_rate);
fprintf('Average length using Huffman on %d-blocks (per sym): %.4f bits/symbol\n', k, avg_len_per_symbol_block);
fprintf('-----------------------------------------------\n');
fprintf('Observation: Huffman(block) and H_k/k should\n');
fprintf('approach H(A) as k increases, confirming\n');
fprintf('better coding efficiency with larger blocks.\n');

end
