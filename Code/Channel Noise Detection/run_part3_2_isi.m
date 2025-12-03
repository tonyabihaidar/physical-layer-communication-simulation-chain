% RUN_PART3_2_ISI
% Part 3.2: Linear Gaussian ISI channel with QPSK and Viterbi MLSE,
% plus AWGN (no-ISI) baseline and theoretical curve.

clear; clc; close all;

rng(2);                            % different seed from Part 3.1 if you like

N_info_bits = 2e5;                % info bits (adjust if runtime too long)
EbN0dB       = 0:2:14;            % Eb/N0 points in dB

% --- ISI channel with Viterbi MLSE (uncoded QPSK) ---
ber_isi_mlse = simulate_qpsk_isi_viterbi(N_info_bits, EbN0dB);

% --- AWGN (no-ISI) QPSK baseline (simulation + theory) ---
[ber_awgn_sim, ber_awgn_theory] = simulate_qpsk_awgn_ebn0(N_info_bits, EbN0dB);

% --- Plot ---
figure;
semilogy(EbN0dB, ber_awgn_theory, 'b-', 'LineWidth', 1.5); hold on;
semilogy(EbN0dB, ber_awgn_sim, 'c--o', 'LineWidth', 1.2);
semilogy(EbN0dB, ber_isi_mlse, 'r-s', 'LineWidth', 1.5);

grid on;
xlabel('E_b/N_0 [dB]');
ylabel('Bit Error Rate (BER)');
legend('AWGN QPSK (theory)', ...
       'AWGN QPSK (sim)', ...
       'ISI + Viterbi MLSE', ...
       'Location','southwest');
title('Part 3.2: QPSK over ISI Channel with Viterbi Equalization');
