%% full_chain.m
%   source x  -> quantizer -> lossless coder -> decoder -> de-quantization
% Measures for each (quantizer, coder):
%   total bits, bits/symbol, latency, SNR, etc.
%
% Also passes the compressed bitstream through a simple
% BPSK AWGN channel, both uncoded and with repetition coding (L=3),
% and stores the resulting channel BER per scheme.

clear; clc;
addpath(genpath('Quantization'));
addpath(genpath('Coding'));
addpath(genpath('Sampling'));
addpath(genpath('Modulation'));
addpath(genpath('Channel Coding'));
%% ---------------- USER PARAMETERS ----------------
Nsymbols = 5e4;     % number of source samples 
fs       = 1;       % "sampling frequency"

% ----- define source (discrete-time) -----
% Here we just generate a zero-mean Gaussian source.
x = randn(Nsymbols,1);
x = x(:);           % column vector
N = numel(x);
signalPower = mean(x.^2);

% ----- Quantizers to test -----
quantizers = { ...
    struct('name','two_level', 'M',2,  'id','Q2'), ...
    struct('name','uniform',   'M',8,  'id','QU8'), ...
    struct('name','uniform',   'M',16, 'id','QU16'), ...
    struct('name','lloydmax',  'M',8,  'id','QL8'), ...
    struct('name','lloydmax',  'M',16, 'id','QL16') ...
};

% ----- Coding strategies -----
% Both have block length = 1 symbol.
coders = { ...
    struct('type','fixed', 'id','CFIX'), ...
    struct('type','huff',  'id','CHUF') ...
};

% ----- Channel parameters for the "physical layer" block -----
% We chose them
EbN0dB_chan = 6;    % example Eb/N0 for channel test (can change)
L_rep       = 3;    % repetition factor for channel coding

results = struct([]);
rIdx    = 0;

%% ---------------- MAIN LOOPS ----------------
for q = 1:numel(quantizers)

    qcfg = quantizers{q};

    % === Quantization ====================================================
    switch qcfg.name
        case 'two_level'
            [xq, thr, lvl, mse_q] = two_level_quant(x);
            M = numel(lvl);

        case 'uniform'
            [xq, thr, lvl, mse_q] = uniform_quant(x, qcfg.M);
            M = numel(lvl);

        case 'lloydmax'
            [xq, thr, lvl, mse_q] = lloydmax_quant(x, qcfg.M);
            M = numel(lvl);

        otherwise
            error('Unknown quantizer type: %s', qcfg.name);
    end

    % SNR due to quantization alone
    SNRdB_q = 10*log10(signalPower / mse_q);

    % Make sure shapes are consistent
    xq = xq(:);

    %% -------- For each coder on this quantized sequence ---------------
    for c = 1:numel(coders)

        ccfg = coders{c};
        rIdx = rIdx + 1;

        schemeName = sprintf('%s_%s', qcfg.id, ccfg.id);

        switch ccfg.type

            case 'fixed'
                % ---------- Fixed-length coding -------------------------
                % Encode (returns bitstring + code->level map)
                tic;
                [bitstr, code2level] = fixedlengthcode_encode(schemeName, xq, lvl);
                encTime = toc;

                % Decode back to levels (ideal, no channel errors)
                tic;
                xq_hat = fixedlengthcode_decode(bitstr, code2level);
                decTime = toc;

                latencySym = 1;          % effectively symbol-by-symbol
                modelTime  = 0;          % no separate modeling stage

            case 'huff'
                % ---------- Symbol-wise Huffman coding ------------------
                % Step 1: build codes + statistics (this also learns p)
                tic;
                [codes, avg_len] = huffman_encode(schemeName, xq, lvl);
                modelTime = toc;         % includes probability estimation & code design

                % Step 2: encode the whole sequence using those codes
                tic;
                bitstr = encode_sequence_huffman(xq, lvl, codes);
                encTime = toc;

                % Step 3: decode back to quantized levels (ideal channel)
                tic;
                xq_hat = huffman_decode(bitstr, codes, lvl);
                decTime = toc;

                latencySym = N;         % we assumed empirical_model uses full sequence

            otherwise
                error('Unknown coder type: %s', ccfg.type);
        end

        % Ensure xq_hat has same shape as x
        xq_hat = xq_hat(:);
        if numel(xq_hat) ~= N
            error('Length mismatch after decoding for scheme %s.', schemeName);
        end

        % === Total distortion and SNR (source + quantizer + lossless) ==
        mse_total   = mean((x - xq_hat).^2);
        SNRdB_total = 10*log10(signalPower / mse_total);

        % === Bit statistics (source coding) =============================
        totalBits   = numel(bitstr);
        bitsPerSym  = totalBits / N;
        latencySec  = latencySym / fs;

        %% === CHANNEL (BPSK over AWGN) ==================================
        % Treat the compressed bitstream "bitstr" as the source for the
        % digital channel. We:
        %   - convert it to {0,1} bits,
        %   - send uncoded BPSK over AWGN,
        %   - send repetition-coded BPSK over AWGN,
        %   - measure BER in both cases.

        [ber_unc_chan, ber_cod_chan] = channel_awgn_bpsk_from_bitstr( ...
            bitstr, EbN0dB_chan, L_rep);

        %% === Store results =============================================
        results(rIdx).quantizerName = qcfg.name;
        results(rIdx).quantizerID   = qcfg.id;
        results(rIdx).M             = M;
        results(rIdx).coderType     = ccfg.type;
        results(rIdx).coderID       = ccfg.id;

        results(rIdx).Nsymbols      = N;
        results(rIdx).signalPower   = signalPower;
        results(rIdx).mse_q         = mse_q;
        results(rIdx).SNRdB_q       = SNRdB_q;
        results(rIdx).mse_total     = mse_total;
        results(rIdx).SNRdB_total   = SNRdB_total;

        results(rIdx).totalBits     = totalBits;
        results(rIdx).bitsPerSym    = bitsPerSym;
        results(rIdx).latencySym    = latencySym;
        results(rIdx).latencySec    = latencySec;
        results(rIdx).modelTime     = modelTime;
        results(rIdx).encTime       = encTime;
        results(rIdx).decTime       = decTime;

        % Channel-related metrics (full chain from source coding outward)
        results(rIdx).EbN0dB_chan      = EbN0dB_chan;
        results(rIdx).L_rep_chan       = L_rep;
        results(rIdx).ber_uncoded_chan = ber_unc_chan;
        results(rIdx).ber_coded_chan   = ber_cod_chan;
    end
end

%% Convert to table and save
T = struct2table(results);
disp(T);

save('full_chain_results.mat', 'results', 'T');
fprintf('\nSaved results to full_chain_results.mat\n');

%% -------- Local helper: encode sequence with Huffman codes ----------
function bitstr = encode_sequence_huffman(xq, lvl, codes)
    % Map each value in xq to the corresponding Huffman code.

    % Make sure column
    xq  = xq(:);
    lvl = lvl(:);

    % Find, for each xq(k), which index in lvl it corresponds to
    [tf, idx] = ismember(xq, lvl);

    % If any value not found in lvl, throw an error
    if ~all(tf)
        error('Some quantized samples are not present in lvl.');
    end

    % codes should be a cell array where codes{i} is a char row like '0101'
    bitCells = codes(idx);      % cell array of code strings
    bitstr   = [bitCells{:}];   % concatenate into a single char row
end


%% -------- Local helper: channel over BPSK/AWGN ----------------------
function [ber_unc, ber_cod] = channel_awgn_bpsk_from_bitstr(bitstr, EbN0dB, L)
    % Convert char '0'/'1' bitstring to numeric bits {0,1}
    b = (bitstr(:) == '1');   % column vector of 0/1

    % Assume Eb = 1 for BPSK => Es = Eb
    EbN0_lin = 10^(EbN0dB/10);
    Eb = 1;
    N0 = Eb / EbN0_lin;
    % For real BPSK, noise variance per real dimension is N0/2
    sigma2 = N0/2;

    % ---------- Uncoded BPSK ----------
    x_unc = bpsk_mod(b, 1);       % A = 1
    n_unc = sqrt(sigma2) * randn(size(x_unc));
    y_unc = x_unc + n_unc;
    b_hat_unc = bpsk_demod(y_unc);
    ber_unc = mean(b ~= b_hat_unc);

    % ---------- Repetition-coded BPSK ----------
    b_enc = rep_encode(b, L);
    x_cod = bpsk_mod(b_enc, 1);
    n_cod = sqrt(sigma2) * randn(size(x_cod));
    y_cod = x_cod + n_cod;
    b_hat_cod_unc = bpsk_demod(y_cod);
    b_hat_cod      = rep_decode(b_hat_cod_unc, L);

    % To be safe, align lengths
    N_common = min(length(b), length(b_hat_cod));
    ber_cod  = mean(b(1:N_common) ~= b_hat_cod(1:N_common));
end
