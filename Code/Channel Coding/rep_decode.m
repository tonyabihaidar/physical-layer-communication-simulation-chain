function b_dec = rep_decode(b_rx, L)
% REP_DECODE   Majority vote decoder for repetition code
%
% Inputs:
%   b_rx : received bits AFTER demod (0/1)
%   L    : repetition factor
%
% Output:
%   b_dec : decoded bitstream (0/1)

    b_rx = b_rx(:);

    % Truncate if last block is incomplete
    N = floor(length(b_rx) / L) * L;
    b_rx = b_rx(1:N);

    % Reshape to L rows: each column is a repeated block
    blocks = reshape(b_rx, L, []);

    % Majority vote on each block
    ones_count = sum(blocks, 1);
    b_dec = ones_count > (L/2);
    b_dec = b_dec(:);

end
