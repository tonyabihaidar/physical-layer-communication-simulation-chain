function [codes, avg_len] = huffman_encode(name, xq, lvl)
%HUFFMAN_CODING  Part 3.2(c) – Huffman coding based on empirical PMF

% Outputs:
%   codes      : cell array of Huffman codewords for each level
%   avg_len    : average codeword length (bits per symbol)

    %--- Step 1: Map levels to indices (1..M)
    [tf, ~] = ismember(xq, lvl);
    if ~all(tf)
        error('Some samples in xq_levels are not found in lvl.');
    end
    M = numel(lvl);

    %--- Step 2: Compute empirical probabilities
    [p,H]=empirical_model(name,xq,lvl);

    %--- Step 3: Build Huffman dictionary using MATLAB built-in function
    dict = huffmandict(num2cell(1:M), p);  % {symbol, code}

    % Extract codes as cell array of strings
    codes = cell(1, M);
    for i = 1:M
        codes{i} = sprintf('%d', dict{i,2});
    end

    %--- Step 4: Compute average code length
    Lvec    = cellfun(@length, codes(:));   % Mx1
    avg_len = sum(p(:) .* Lvec);            % scalar


    %--- Step 5: Display results
    fprintf('\n--- %s: Huffman Coding ---\n', name);
    fprintf('Level\tp(level)\tCode\n');
    fprintf('--------------------------------\n');
    for i = 1:M
        fprintf('%6.2f\t%.4f\t\t%s\n', lvl(i), p(i), codes{i});
    end
    fprintf('--------------------------------\n');
    fprintf('Average code length = %.4f bits/symbol\n', avg_len);
end
