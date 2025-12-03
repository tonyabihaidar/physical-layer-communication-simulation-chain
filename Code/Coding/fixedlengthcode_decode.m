function xq = fixedlengthcode_decode(bitstr, code2level)
% Decode fixed-length bitstring back to *levels* using code→level map.
% bitstr    : char array of '0'/'1'
% code2level: containers.Map('char' -> level) from flc_make_code_map

    ks = keys(code2level);
    if isempty(ks), error('Empty code2level map.'); end
    L = length(ks{1});                 % number of bits per symbol
    N = floor(numel(bitstr)/L);        % number of symbols

    xq = zeros(1, N);
    pos = 1;
    for n = 1:N
        code = bitstr(pos:pos+L-1);
        if ~isKey(code2level, code)
            error('Unknown code "%s" at symbol %d.', code, n);
        end
        xq(n) = code2level(code);
        pos = pos + L;
    end
end
