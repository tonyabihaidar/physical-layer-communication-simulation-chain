% RUN_PART3_1_AWGN
% Full Part 3.1: AWGN channel, coded and uncoded, BPSK and QPSK.

clear; clc; close all;
addpath(genpath('../Modulation'));
addpath(genpath('../Channel Coding'));

rng(1);                             % reproducibility

N_info_bits = 1e6;                  % number of information bits (before coding)
SNRdB       = 0:2:12;               % SNR points in dB: SNR = 10 log10(Es)
L_rep       = 3;                    % repetition factor for channel coding

% ---------- BPSK ----------
[ber_bpsk_unc, ber_bpsk_cod] = simulate_bpsk_awgn_coded(N_info_bits, SNRdB, L_rep);

% ---------- QPSK ----------
[ber_qpsk_unc, ber_qpsk_cod] = simulate_qpsk_awgn_coded(N_info_bits, SNRdB, L_rep);

% ---------- Plot ----------
figure;
semilogy(SNRdB, ber_bpsk_unc, '-o', 'LineWidth', 1.5); hold on;
semilogy(SNRdB, ber_bpsk_cod, '-o', 'LineWidth', 1.5, 'MarkerFaceColor','auto');
semilogy(SNRdB, ber_qpsk_unc, '-s', 'LineWidth', 1.5);
semilogy(SNRdB, ber_qpsk_cod, '-s', 'LineWidth', 1.5, 'MarkerFaceColor','auto');
grid on;
xlabel('SNR = 10 log_{10}(E_s) [dB]');
ylabel('Bit Error Rate (BER)');
legend('BPSK uncoded', ...
       ['BPSK coded (L = ' num2str(L_rep) ')'], ...
       'QPSK uncoded', ...
       ['QPSK coded (L = ' num2str(L_rep) ')'], ...
       'Location','southwest');
title('Part 3.1: Coded vs Uncoded BER over AWGN');
