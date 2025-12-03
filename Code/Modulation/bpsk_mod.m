function x = bpsk_mod(b, A)
% BPSK_MOD: map bits {0,1} to BPSK symbols {+A, -A}

    b = b(:);             % ensure column
    x = A * (2*b - 1);    % 0 -> -A, 1 -> +A

end
