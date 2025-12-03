%% part3_2_3_ber_sim.m
% Part 3.2.3 – BER Estimation Procedure (Simulation)
% QPSK over:
%   - ISI channel: y_i = 2 x_i + x_{i-1} + n_i, detected by Viterbi (MLSE)
%   - AWGN channel: y_i = x_i + n_i
% And compare to the theoretical AWGN QPSK BER: Q(sqrt(2 Eb/N0))

clear; clc; close all;
addpath(genpath('../Modulation'));

rng(3);                        % seed for reproducibility

% ---- Simulation parameters ----
Nbits   = 2e5;                 % number of info bits (must be even)
EbN0dB  = 0:2:14;              % Eb/N0 points in dB

% Ensure even number of bits for QPSK
if mod(Nbits,2) ~= 0
    Nbits = Nbits - 1;
end

% Preallocate
ber_isi_mlse   = zeros(size(EbN0dB));
ber_awgn_sim   = zeros(size(EbN0dB));
ber_awgn_theor = zeros(size(EbN0dB));

% QPSK parameters: Es = 1 => A = 1
Es = 1;
A  = sqrt(Es);                 % = 1

% ---- (a) Generate bits (once) ----
b_info = randi([0 1], Nbits, 1);

fprintf('=== Part 3.2.3: BER over ISI (Viterbi) and AWGN (QPSK) ===\n');

for k = 1:length(EbN0dB)

    % ---- Eb/N0 and noise variance ----
    EbN0_lin = 10^(EbN0dB(k)/10);
    Eb       = 1;               % as specified in the project
    N0       = Eb / EbN0_lin;   % N0 = 1 / (Eb/N0)
    sigma2   = N0;              % variance per complex symbol

    %% ========== ISI CHANNEL + VITERBI (MLSE) ==========
    % (a) Map bits to QPSK symbols x in S, Es = 1
    x = qpsk_mod(b_info, A);    % complex QPSK, Es = 1

    Nsym = length(x);
    y_isi = zeros(Nsym,1);

    % (b) ISI channel: y_i = 2 x_i + x_{i-1} + n_i
    n_isi = sqrt(sigma2/2) .* (randn(Nsym,1) + 1j*randn(Nsym,1));
    x_prev = 0;                 % initial previous symbol

    for i = 1:Nsym
        y_isi(i) = 2*x(i) + x_prev + n_isi(i);
        x_prev   = x(i);
    end

    % (c) MLSE with Viterbi equalizer
    x_hat_isi = viterbi_equalizer(y_isi);   % estimated QPSK symbols
    b_hat_isi = qpsk_demod(x_hat_isi);      % back to bits

    ber_isi_mlse(k) = mean(b_info ~= b_hat_isi);

    %% ========== AWGN BASELINE (NO ISI) ==========
    % (d) y_i = x_i + n_i over AWGN only
    n_awgn = sqrt(sigma2/2) .* (randn(Nsym,1) + 1j*randn(Nsym,1));
    y_awgn = x + n_awgn;

    b_hat_awgn = qpsk_demod(y_awgn);
    ber_awgn_sim(k) = mean(b_info ~= b_hat_awgn);

    % Theoretical AWGN QPSK BER: Q( sqrt(2 Eb/N0) )
    % Q(x) = 0.5 * erfc(x / sqrt(2)) => Q(sqrt(2*gamma)) = 0.5 * erfc(sqrt(gamma))
    ber_awgn_theor(k) = 0.5 * erfc( sqrt(EbN0_lin) );

    fprintf('Eb/N0 = %2d dB:  BER_ISI+MLSE = %.3e,  AWGN_sim = %.3e,  AWGN_theor = %.3e\n', ...
        EbN0dB(k), ber_isi_mlse(k), ber_awgn_sim(k), ber_awgn_theor(k));
end

%% ========== PLOT RESULTS ==========
figure;
semilogy(EbN0dB, ber_awgn_theor, 'k-',  'LineWidth', 1.5); hold on;
semilogy(EbN0dB, ber_awgn_sim,   'ko--','LineWidth', 1.2);
semilogy(EbN0dB, ber_isi_mlse,   'rs-','LineWidth', 1.5);

grid on;
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('AWGN QPSK (theory)', ...
       'AWGN QPSK (sim)', ...
       'ISI + Viterbi MLSE', ...
       'Location','southwest');
title('Part 3.2.3: BER for QPSK over ISI (MLSE) and AWGN');
