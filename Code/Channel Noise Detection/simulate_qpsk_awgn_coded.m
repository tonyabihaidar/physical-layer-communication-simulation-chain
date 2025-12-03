function [ber_unc, ber_cod] = simulate_qpsk_awgn_coded(N_info_bits, SNRdB, L)
% SIMULATE_QPSK_AWGN_CODED
%   QPSK over AWGN, uncoded and repetition-coded.
%
% Inputs:
%   N_info_bits : number of information bits (before coding)
%   SNRdB       : vector of SNR values (SNR = 10 log10(Es))
%   L           : repetition factor
%
% Outputs:
%   ber_unc : uncoded BER vs SNR
%   ber_cod : coded BER vs SNR

    % We need an even number of bits for QPSK (pairs).
    % Also want N_info_bits * L to be even for coded case.
    % We'll adjust by at most 1 bit if needed (negligible).

    N_info_bits_adj = N_info_bits;

    % Make N_info_bits even so uncoded QPSK has even number of bits
    if mod(N_info_bits_adj, 2) ~= 0
        N_info_bits_adj = N_info_bits_adj - 1;
    end

    % Ensure N_info_bits_adj * L is even
    if mod(N_info_bits_adj * L, 2) ~= 0
        N_info_bits_adj = N_info_bits_adj - 1;
    end

    % Generate final information bits used in both paths
    b_info = randi([0 1], N_info_bits_adj, 1);

    ber_unc = zeros(size(SNRdB));
    ber_cod = zeros(size(SNRdB));

    fprintf('=== QPSK over AWGN (uncoded + coded, L = %d) ===\n', L);

    for k = 1:length(SNRdB)
        Es = 10^(SNRdB(k)/10);   % symbol energy
        A  = sqrt(Es);           % amplitude

        % ---------- Uncoded path ----------
        b_tx_unc = b_info;       % already even length
        x_unc    = qpsk_mod(b_tx_unc, A);  % complex QPSK

        Nsym_unc = length(x_unc);
        n_unc = sqrt(0.5) * (randn(Nsym_unc,1) + 1j*randn(Nsym_unc,1));
        y_unc = x_unc + n_unc;

        b_hat_unc = qpsk_demod(y_unc);
        ber_unc(k) = mean(b_info ~= b_hat_unc);

        % ---------- Coded path (repetition L) ----------
        b_tx_cod = rep_encode(b_info, L);       % length is N_info_bits_adj * L

        % Just in case of any odd length, ensure even length for QPSK
        if mod(length(b_tx_cod),2) ~= 0
            b_tx_cod = b_tx_cod(1:end-1);
        end

        x_cod = qpsk_mod(b_tx_cod, A);

        Nsym_cod = length(x_cod);
        n_cod = sqrt(0.5) * (randn(Nsym_cod,1) + 1j*randn(Nsym_cod,1));
        y_cod = x_cod + n_cod;

        b_hat_cod_rx = qpsk_demod(y_cod);       % raw decoded coded bits

        % rep_decode expects a length multiple of L
        N_rx = floor(length(b_hat_cod_rx)/L) * L;
        b_hat_cod_rx = b_hat_cod_rx(1:N_rx);
        b_hat_cod     = rep_decode(b_hat_cod_rx, L);

        % b_hat_cod length should match b_info (possibly minus 0–1 bit);
        N_common = min(length(b_info), length(b_hat_cod));

        ber_cod(k) = mean(b_info(1:N_common) ~= b_hat_cod(1:N_common));

        fprintf('SNR = %2d dB:  BER_unc = %.3e,  BER_cod = %.3e\n', ...
            SNRdB(k), ber_unc(k), ber_cod(k));
    end
end
