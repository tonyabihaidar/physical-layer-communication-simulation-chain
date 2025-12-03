%% analyze_results.m
% Uses results saved by full_chain.m
% (i)   Builds compact table with total bits, bits/symbol, latency + NEW channel BER fields.
% (j)   Finds useful SNR-matched comparisons.
% (k)   NEW: Adds channel BER analysis (uncoded vs coded).

clear; clc;

load('full_chain_results.mat', 'results', 'T');

%% ---------------------------------------------------------------
% 1) Compact summary table (with channel BER added)
%% ---------------------------------------------------------------
fprintf('===== End-to-end metrics for each (quantizer, coder) =====\n');

summaryTable = T(:, { ...
    'quantizerID', 'coderID', 'M', ...
    'SNRdB_total', 'totalBits', 'bitsPerSym', ...
    'latencySym', 'latencySec', ...
    'EbN0dB_chan', 'ber_uncoded_chan', 'ber_coded_chan'});

disp(summaryTable);

%% ---------------------------------------------------------------
% 2) Bits/symbol per (quantizer, coder)
%% ---------------------------------------------------------------
fprintf('\n===== Average bits/symbol per quantizer & coder =====\n');

uQ = unique(T.quantizerID, 'stable');
uC = unique(T.coderID,     'stable');

avgBits = nan(numel(uQ), numel(uC));

for i = 1:numel(uQ)
    for j = 1:numel(uC)
        mask = strcmp(T.quantizerID, uQ{i}) & strcmp(T.coderID, uC{j});
        if any(mask)
            avgBits(i,j) = mean(T.bitsPerSym(mask));
        end
    end
end

bitsTable = array2table(avgBits, 'RowNames', uQ, 'VariableNames', uC);
disp(bitsTable);

%% ---------------------------------------------------------------
% 3) (j) Required comparisons at equal SNR
%% ---------------------------------------------------------------
tolSNR = 0.1;   % dB tolerance for "equal SNR"

%% 3A) Same quantizer, different coders
fprintf('\n===== Same quantizer, different coders (identical SNR) =====\n');

rowsA = {};
SNR_A = [];
bitsA = [];

for i = 1:height(T)
    for j = i+1:height(T)
        if strcmp(T.quantizerID{i}, T.quantizerID{j}) && ...
           ~strcmp(T.coderID{i}, T.coderID{j}) && ...
           abs(T.SNRdB_total(i) - T.SNRdB_total(j)) <= tolSNR

            rowsA{end+1,1} = sprintf('Q=%s: %s vs %s', ...
                T.quantizerID{i}, T.coderID{i}, T.coderID{j}); %#ok<AGROW>

            SNR_A(end+1,1)  = T.SNRdB_total(i); %#ok<AGROW>
            bitsA(end+1,1)  = T.bitsPerSym(i); %#ok<AGROW>
            bitsA(end,2)    = T.bitsPerSym(j); %#ok<AGROW>
        end
    end
end

if isempty(rowsA)
    fprintf('No such pairs found.\n');
else
    TA = table(rowsA, SNR_A, bitsA(:,1), bitsA(:,2), ...
        'VariableNames', {'Pair','SNRdB','bitsPerSym_1','bitsPerSym_2'});
    disp(TA);
end

%% 3B) Same coder, different quantizers with nearly equal SNR
fprintf('\n===== Same coder, different quantizers (≈ equal SNR) =====\n');

rowsB = {};
SNR1 = []; SNR2 = [];
bits1 = []; bits2 = [];

for i = 1:height(T)
    for j = i+1:height(T)
        if strcmp(T.coderID{i}, T.coderID{j}) && ...
           ~strcmp(T.quantizerID{i}, T.quantizerID{j}) && ...
           abs(T.SNRdB_total(i) - T.SNRdB_total(j)) <= tolSNR

            rowsB{end+1,1} = sprintf('%s_%s  vs  %s_%s', ...
                T.quantizerID{i}, T.coderID{i}, ...
                T.quantizerID{j}, T.coderID{j}); %#ok<AGROW>

            SNR1(end+1,1)  = T.SNRdB_total(i); %#ok<AGROW>
            SNR2(end+1,1)  = T.SNRdB_total(j); %#ok<AGROW>
            bits1(end+1,1) = T.bitsPerSym(i);  %#ok<AGROW>
            bits2(end+1,1) = T.bitsPerSym(j);  %#ok<AGROW>
        end
    end
end

if isempty(rowsB)
    fprintf('No SNR-matched quantizer pairs with tol = %.2f dB.\n', tolSNR);
else
    TB = table(rowsB, SNR1, SNR2, bits1, bits2, ...
        'VariableNames', {'Pair','SNR1_dB','SNR2_dB','bitsPerSym1','bitsPerSym2'});
    disp(TB);
end

%% ---------------------------------------------------------------
% 4) NEW: Channel BER summary table
%% ---------------------------------------------------------------
fprintf('\n===== Channel BER summary (uncoded vs coded) =====\n');

chanTable = T(:, { ...
    'quantizerID', 'coderID', ...
    'EbN0dB_chan', 'ber_uncoded_chan', 'ber_coded_chan'});

disp(chanTable);

%% ---------------------------------------------------------------
% 5) NEW: Channel BER improvement (coding gain table)
%% ---------------------------------------------------------------
fprintf('\n===== Channel Coding Gain (BER_uncoded / BER_coded) =====\n');

gain = T.ber_uncoded_chan ./ T.ber_coded_chan;
gainTable = table(T.quantizerID, T.coderID, T.EbN0dB_chan, ...
    T.ber_uncoded_chan, T.ber_coded_chan, gain, ...
    'VariableNames', {'Quantizer','Coder','EbN0dB','BER_uncoded','BER_coded','Gain'});

disp(gainTable);

%% End of script
fprintf('\nAnalysis completed.\n');
