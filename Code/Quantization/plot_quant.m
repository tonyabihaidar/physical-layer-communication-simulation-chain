function plot_quant(t, x, xq, thr, lvl, title_text)
%PLOT_QUANT  Generic plotting utility for quantization results.
%
%   plot_quant(t, x, xq, thr, lvl, title_text)
%   - t : time or sample indices (vector same length as x)
%   - x : original signal
%   - xq: quantized signal (from quan)
%   - thr: scalar or vector of thresholds
%   - lvl: vector of quantization levels
%   - title_text: title for the plot

    figure;
    plot(t, x, 'LineWidth', 1.3); hold on;
    stairs(t, xq, 'LineWidth', 1.1);

    % Plot thresholds as dashed horizontal lines
    if ~isempty(thr)
        for i = 1:length(thr)
            yline(thr(i), '--', sprintf('thr%d', i), 'LabelHorizontalAlignment', 'left');
        end
    end

    grid on;
    legend('Original', 'Quantized', 'Location', 'best');
    xlabel('Samples or time'); ylabel('Amplitude');
    title(title_text);

    % Show levels numerically in console
    disp('Quantization levels:');
    disp(lvl);
end
