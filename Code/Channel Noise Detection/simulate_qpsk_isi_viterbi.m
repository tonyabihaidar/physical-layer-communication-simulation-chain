function ber_isi = simulate_qpsk_isi_viterbi(N_info_bits, EbN0dB)
% SIMULATE_QPSK_ISI_VITERBI
%   Uncoded QPSK over ISI channel y_i = 2 x_i + x_{i-1} + n_i
%   Detection via Viterbi MLSE, BER vs Eb/N0.
%
% Inputs:
%   N_info_bits : number of information bits (must be even)
%   EbN0dB      : vector of Eb/N0 values in dB
%
% Output:
%   ber_isi : BER vs Eb/N0 for ISI + Viterbi MLSE

    % Ensure even number of bits (2 bits per QPSK symbol)
    if mod(N_info_bits, 2) ~= 0
        N_info_bits = N_info_bits - 1;
    end

    % Generate information bits
    b_info = randi([0 1], N_info_bits, 1);

    ber_isi = zeros(size(EbN0dB));

    % --- QPSK parameters ---
    Es = 1;           % symbol energy
    A  = sqrt(Es);    % amplitude (== 1 here)

    fprintf('=== ISI channel + Viterbi MLSE (QPSK) ===\n');

    for k = 1:length(EbN0dB)

        % Eb/N0 (linear)
        EbN0_lin = 10^(EbN0dB(k)/10);

        % Using Eb = 1 -> N0 = 1 / (Eb/N0)
        N0      = 1 / EbN0_lin;
        sigma2  = N0;                 % per complex symbol

        % ---- Modulate bits to QPSK symbols ----
        x = qpsk_mod(b_info, A);      % complex QPSK symbols, Es = 1

        Nsym = length(x);
        y    = zeros(Nsym,1);

        % ---- ISI channel: y_i = 2 x_i + x_{i-1} + n_i ----
        x_prev = 0;   % initial previous symbol (can be 0; effect negligible for long N)

        % Generate complex AWGN with variance sigma2 per symbol
        n = sqrt(sigma2/2) .* (randn(Nsym,1) + 1j*randn(Nsym,1));

        for i = 1:Nsym
            y(i)   = 2*x(i) + x_prev + n(i);
            x_prev = x(i);
        end

        % ---- Viterbi equalizer (MLSE) ----
        x_hat = viterbi_equalizer(y);

        % ---- Demap symbols back to bits ----
        b_hat = qpsk_demod(x_hat);

        % BER
        ber_isi(k) = mean(b_info ~= b_hat);

        fprintf('Eb/N0 = %2d dB:  BER_ISI+MLSE = %.3e\n', EbN0dB(k), ber_isi(k));
    end

end
