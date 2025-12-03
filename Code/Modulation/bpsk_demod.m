function b_hat = bpsk_demod(y)
% BPSK_DEMOD: threshold at zero

    b_hat = y >= 0;

end
