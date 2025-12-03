function x = qpsk_mod(b, A)
% QPSK_MOD: Gray-mapped QPSK

    b = b(:);
    if mod(length(b),2) ~= 0
        error('QPSK requires even number of bits.');
    end

    b2 = reshape(b,2,[]).';  % group bits into pairs
    Nsym = size(b2,1);
    x = zeros(Nsym,1);
    s = A/sqrt(2);

    for i = 1:Nsym
        pair = b2(i,:);
        if     all(pair == [1 1])
            x(i) = s*( 1 + 1j);
        elseif all(pair == [1 0])
            x(i) = s*(-1 + 1j);
        elseif all(pair == [0 0])
            x(i) = s*(-1 - 1j);
        elseif all(pair == [0 1])
            x(i) = s*( 1 - 1j);
        end
    end
end
