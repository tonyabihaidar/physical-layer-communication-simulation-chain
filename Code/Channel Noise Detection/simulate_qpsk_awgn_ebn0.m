function [ber_awgn, ber_theory] = simulate_qpsk_awgn_ebn0(N_info_bits, EbN0dB)
% SIMULATE_QPSK_AWGN_EBN0
%   Uncoded QPSK over AWGN (no ISI), BER vs Eb/N0.
%   Also returns the theoretical BER curve Q(sqrt(2 Eb/N0)).
%
% Inputs:
%   N_info_bits : number of information bits (must be even)
%   EbN0dB      : vector of Eb/N0 values in dB
%
% Outputs:
%   ber_awgn   : simulated BER vs Eb/N0
%   ber_theory : theoretical BER vs Eb/N0

    % Ensure even bits
    if mod(N_info_bits, 2) ~= 0
        N_info_bits = N_info_bits - 1;
    end

    % Info bits
    b_info = randi([0 1], N_info_bits, 1);

    ber_awgn   = zeros(size(EbN0dB));
    ber_theory = zeros(size(EbN0dB));

    % QPSK parameters
    Es = 1;           % symbol energy
    A  = sqrt(Es);    % amplitude = 1

    fprintf('=== AWGN QPSK baseline (no ISI) ===\n');

    for k = 1:length(EbN0dB)

        EbN0_lin = 10^(EbN0dB(k)/10);
        N0       = 1 / EbN0_lin;        % with Eb = 1
        sigma2   = N0;                  % per complex symbol

        % Modulation
        x = qpsk_mod(b_info, A);

        Nsym = length(x);
        n    = sqrt(sigma2/2) .* (randn(Nsym,1) + 1j*randn(Nsym,1));
        y    = x + n;

        % Demodulation
        b_hat = qpsk_demod(y);
        ber_awgn(k) = mean(b_info ~= b_hat);

        % Theoretical BER: Q(sqrt(2 Eb/N0))
        % => Q(sqrt(2*gamma)) = 0.5 * erfc(sqrt(gamma))
        ber_theory(k) = 0.5 * erfc(sqrt(EbN0_lin));

        fprintf('Eb/N0 = %2d dB:  BER_AWGN_sim = %.3e,  BER_theory = %.3e\n', ...
            EbN0dB(k), ber_awgn(k), ber_theory(k));
    end

end
