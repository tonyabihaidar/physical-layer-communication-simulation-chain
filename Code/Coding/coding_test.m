%% Part 3 Compare Fixed-Length, Huffman, and Block Coding


clear; clc;

%% ----- Inputs: replace with your actual quantizer outputs -----
name = 'Quantizer A';
seed = 442;

% Example levels and a toy stream (LEVELS, not indices):
lvl = [-3.6 -1.2 1 3.4 5.6];
xq  = [-3.6 -1.2 -1.2 1 5.6 3.4 3.4 -3.6 -1.2 1 3.4 3.4 5.6 -3.6 1 3.4 -1.2 -3.6 -1.2 -1.2 1 5.6 3.4 3.4 -3.6 -1.2 1 3.4 3.4 5.6 -3.6 1 3.4 -1.2 -3.6 -1.2 -1.2 1 5.6 3.4 3.4 -3.6 -1.2 1 3.4 3.4 5.6 -3.6 1 3.4 -1.2 -3.6 -1.2 -1.2 1 5.6 3.4 3.4 -3.6 -1.2 1 3.4 3.4 5.6 -3.6 1 3.4 -1.2];

% (Optional) thresholds just for display; replace if you have them
thr = [-2 -0.5 0.5 2];

%% 0) Part 3.1 display + fix seed
display_quantizer_info(name, xq, thr, lvl, seed);

%% 1) Empirical model (p and H)
[p, H] = empirical_model(name, xq, lvl);

%% 2) Fixed-length coding (Unknown distribution pass)
[bitstr_fl, code2level_fl] = fixedlengthcode_encode([name ' (FLC)'], xq, lvl);
xq_fl = fixedlengthcode_decode(bitstr_fl, code2level_fl);

% Metrics
M     = numel(lvl);
Lfix  = ceil(log2(M));
N     = numel(xq);
bits_fl_total = numel(bitstr_fl);
bps_fl = bits_fl_total / N;

%% 3) Huffman coding (build codes, pack, decode)
[codes, Lavg_huff] = huffman_encode([name ' (Huffman)'], xq, lvl);
bitstr_huff = huffman_pack_bits(xq, lvl, codes);
xq_huff = huffman_decode(bitstr_huff, codes, lvl);

% Metrics
bits_huff_total = numel(bitstr_huff);
bps_huff = bits_huff_total / N;

%% 4) Checks
ok_fl   = isequal(xq, xq_fl);
ok_huff = isequal(xq, xq_huff);

%% 5) Summary printout
fprintf('\n==================== Part 3(d) Summary ====================\n');
fprintf('Entropy H(A)               : %.4f bits/symbol\n', H);
fprintf('Fixed-length L (ceil log2M): %d bits/symbol\n', Lfix);
fprintf('Huffman avg length         : %.4f bits/symbol\n', Lavg_huff);
fprintf('-----------------------------------------------------------\n');
fprintf('Sequence length N          : %d symbols\n', N);
fprintf('Fixed-length total bits    : %d  (%.4f bits/sym)\n', bits_fl_total, bps_fl);
fprintf('Huffman total bits         : %d  (%.4f bits/sym)\n', bits_huff_total, bps_huff);
fprintf('Lossless (Fixed-length)    : %s\n', ternary(ok_fl,'PASS','FAIL'));
fprintf('Lossless (Huffman)         : %s\n', ternary(ok_huff,'PASS','FAIL'));
fprintf('===========================================================\n');

%% 6) First few encoded symbols for each method
Kshow = min(8, N);
L = Lfix;
fprintf('\nFirst %d FLC symbols (L=%d):\n', Kshow, L);
show_first_codes_fixed(xq(1:Kshow), lvl, L);

fprintf('\nFirst %d Huffman symbols:\n', Kshow);
show_first_codes_huff(xq(1:Kshow), lvl, codes);

%% 7) Block coding analysis (parts e, f, g)
fprintf('\n================ Block Coding Analysis ================\n');
for k = [2 3 4]
    if N < k
        fprintf('Skipping k =%d (N<k)\n', k);
        continue;
    end
    % Let the function handle all the detailed prints for part (g)
    fprintf('Concatenation k= %d\n',k);
    block_coding_analysis(name, xq, lvl, k, H, Lavg_huff);
end
fprintf('======================================================\n');

%% -------------------- Local helper functions --------------------
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function bitstr = huffman_pack_bits(xq, lvl, codes)
% Concatenate Huffman codewords (strings) for the LEVEL sequence xq
    [tf, idx] = ismember(xq, lvl);
    if ~all(tf), error('Some samples in xq not found in lvl.'); end
    parts = cell(1, numel(idx));
    for n = 1:numel(idx)
        parts{n} = codes{idx(n)};
    end
    bitstr = strjoin(parts, '');
end

function show_first_codes_fixed(xq_head, lvl, L)
% Print the fixed-length codes for the first few LEVEL symbols
    [tf, idx] = ismember(xq_head, lvl);
    if ~all(tf), error('show_first_codes_fixed: sample not in lvl'); end
    fprintf('Level\tCode\n');
    for i = 1:numel(idx)
        fprintf('%6.2f\t%s\n', xq_head(i), dec2bin(idx(i)-1, L));
    end
end

function show_first_codes_huff(xq_head, lvl, codes)
% Print the Huffman codes for the first few LEVEL symbols
    [tf, idx] = ismember(xq_head, lvl);
    if ~all(tf), error('show_first_codes_huff: sample not in lvl'); end
    fprintf('Level\tCode\n');
    for i = 1:numel(idx)
        fprintf('%6.2f\t%s\n', xq_head(i), codes{idx(i)});
    end
end
