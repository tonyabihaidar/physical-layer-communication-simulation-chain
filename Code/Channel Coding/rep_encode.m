function b_enc = rep_encode(b, L)
% REP_ENCODE   Repetition coding (rate = 1/L)
%
% Inputs:
%   b : input bit column vector (0/1)
%   L : repetition factor
%
% Output:
%   b_enc : encoded bitstream (each bit repeated L times)

    b = b(:);  % ensure column
    b_enc = repelem(b, L);

end
