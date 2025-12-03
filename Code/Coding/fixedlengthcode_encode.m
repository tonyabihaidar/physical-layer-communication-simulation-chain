function [bitstr, code2level] = fixedlengthcode_encode(name, xq, lvl)
    M = length(lvl);
    L = ceil(log2(M));
    N = numel(xq);

    fprintf('\n--- %s: Fixed-Length Coding ---\n', name);
    fprintf('Alphabet size M = %d\n', M);
    fprintf('Bits per symbol L = %d\n', L);

    [tf, idx] = ismember(xq, lvl);
    if ~all(tf)
        error('Some samples are not in lvl. Check quantizer levels vs input.');
    end

    % Build full code -> level map once
    code2level = containers.Map('KeyType','char','ValueType','double');
    for j = 1:M
        code2level(dec2bin(j-1, L)) = lvl(j);
    end

    % Encode sequence
    bitstr = repmat('0', 1, N*L);
    pos = 1;
    for i = 1:N
        b = dec2bin(idx(i)-1, L);
        bitstr(pos:pos+L-1) = b;
        pos = pos + L;
    end
end
