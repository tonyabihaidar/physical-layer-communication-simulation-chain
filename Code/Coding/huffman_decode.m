function xq_recovered = huffman_decode(bitstr, codes, lvl)
%HUFFMAN_DECODE  Decode a Huffman bitstring back to LEVELS.
% Inputs:
%   bitstr : char array of '0'/'1' (concatenated Huffman codes)
%   codes  : 1xM cell array of codewords (strings), aligned with lvl
%   lvl    : 1xM levels (alphabet) in the same order as 'codes'
%
% Output:
%   xq_recovered : 1xN recovered sequence in LEVELS (not indices)

    % Build code -> level map
    M = numel(lvl);
    code2level = containers.Map('KeyType','char','ValueType','double');
    for i = 1:M
        code2level(codes{i}) = lvl(i);
    end

    % Greedy prefix decode
    xq_recovered = [];
    buf = "";
    for k = 1:numel(bitstr)
        buf = buf + bitstr(k);
        key = char(buf);
        if isKey(code2level, key)
            xq_recovered(end+1) = code2level(key);
            buf = ""; % reset buffer
        end
    end

    if strlength(buf) ~= 0
        error('Decoding error: leftover bits "%s" do not match any code.', char(buf));
    end
end
