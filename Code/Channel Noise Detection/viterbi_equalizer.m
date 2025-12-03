function x_hat = viterbi_equalizer(y)
% VITERBI_EQUALIZER
%   ML sequence detector for ISI channel:
%       y_i = 2 s_i + s_{i-1} + n_i
%   with QPSK constellation S = (1/sqrt(2)) * (±1 ± j).
%
% Input:
%   y : received complex sequence (column vector)
%
% Output:
%   x_hat : estimated transmitted QPSK symbols (column vector)

    % ---- QPSK constellation S ----
    S = [ 1+1j; -1+1j; -1-1j; 1-1j ] / sqrt(2);
    Ns = length(S);     % number of states = 4
    N  = length(y);     % time length

    % ---- Path metrics initialization ----
    PM = inf(Ns,1);
    PM(1) = 0;          % start state (arbitrary choice; effect negligible for large N)

    % ---- Storage for traceback ----
    prev_state  = zeros(Ns, N);   % from which state we came
    prev_symbol = zeros(Ns, N);   % which symbol index was chosen

    % ---- Forward recursion ----
    for i = 1:N

        PM_new = inf(Ns,1);

        for sp = 1:Ns            % previous state (symbol index for s_{i-1})
            for sc = 1:Ns        % current symbol index for s_i

                s_prev = S(sp);
                s_curr = S(sc);

                % channel model: y_i ≈ 2 s_i + s_{i-1}
                y_pred = 2*s_curr + s_prev;

                branch_metric = abs( y(i) - y_pred )^2;
                candidate_metric = PM(sp) + branch_metric;

                % keep best path entering state "sc"
                if candidate_metric < PM_new(sc)
                    PM_new(sc)       = candidate_metric;
                    prev_state(sc,i) = sp;
                    prev_symbol(sc,i)= sc;
                end
            end
        end

        PM = PM_new;
    end

    % ---- Traceback ----
    % Start from state with minimum final metric
    [~, state] = min(PM);
    sym_idx_hat = zeros(N,1);

    for i = N:-1:1
        sym_idx_hat(i) = prev_symbol(state, i);
        state          = prev_state(state, i);
        if state == 0
            state = 1;   % safety in case of first step
        end
    end

    % Map indices back to actual QPSK symbols
    x_hat = S(sym_idx_hat);

end
