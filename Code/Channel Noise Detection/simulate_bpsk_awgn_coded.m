function [ber_unc, ber_cod] = simulate_bpsk_awgn_coded(N_info_bits, SNRdB, L)
% SIMULATE_BPSK_AWGN_CODED
%   BPSK over AWGN, uncoded and repetition-coded.
%
% Inputs:
%   N_info_bits : number of information bits (before coding)
%   SNRdB       : vector of SNR values (SNR = 10 log10(Es))
%   L           : repetition factor (e.g., 3)
%
% Outputs:
%   ber_unc : uncoded BER vs SNR
%   ber_cod : coded BER vs SNR

    % Generate information bits (same used for uncoded and coded)
    b_info = randi([0 1], N_info_bits, 1);

    ber_unc = zeros(size(SNRdB));
    ber_cod = zeros(size(SNRdB));

    fprintf('=== BPSK over AWGN (uncoded + coded, L = %d) ===\n', L);

    for k = 1:length(SNRdB)
        Es = 10^(SNRdB(k)/10);   % symbol energy
        A  = sqrt(Es);           % amplitude

        % ----- Uncoded path -----
        b_tx_unc = b_info;
        x_unc    = bpsk_mod(b_tx_unc, A);

        n_unc = sqrt(0.5) * randn(size(x_unc));   % real AWGN, var=0.5
        y_unc = x_unc + n_unc;

        b_hat_unc = bpsk_demod(y_unc);
        ber_unc(k) = mean(b_info ~= b_hat_unc);

        % ----- Coded path (repetition L) -----
        b_tx_cod = rep_encode(b_info, L);         % encode
        x_cod    = bpsk_mod(b_tx_cod, A);

        n_cod = sqrt(0.5) * randn(size(x_cod));
        y_cod = x_cod + n_cod;

        b_hat_cod_rx = bpsk_demod(y_cod);        % hard decisions on coded bits
        b_hat_cod    = rep_decode(b_hat_cod_rx, L);  % majority vote decode

        % BER w.r.t. original information bits
        ber_cod(k) = mean(b_info ~= b_hat_cod);

        fprintf('SNR = %2d dB:  BER_unc = %.3e,  BER_cod = %.3e\n', ...
            SNRdB(k), ber_unc(k), ber_cod(k));
    end
end
