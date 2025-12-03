function [t_sample, x_sample] = sample(t, xt, fs)

    Ts = 1 / fs;
    t_sample = t(1) : Ts : t(end);
    %here we generated a vector in order to sample t which will be further used
    %in order to compute the x_samples
    %using interpolation which is a function in MATLAB
    x_sample = interp1(t, xt, t_sample);

end

function xrcon = reconstruct(t, x_sample, fs)

    Ts = 1 / fs;
    t_sample = t(1) : Ts : t(end); 
    %again, same basis here, we sampled t and that will be used in order
    %to compute xrcon using interpolation, the mode "spline" was used
    %to give the function a smooth curve when connecting the values
    xrcon = interp1(t_sample, x_sample, t, 'spline');
end


%Now to start, we'll define the signal and all of its values
f0 = 1000;              % Signal frequency
fN = 2 * f0;            % Nyquist rate
fs_orig = 100 * fN;     % Very high sampling rate for a "proper"cos function to compare to
T_end = 2 * (1/f0);     %here I made it 2* so we can better the graph

t = 0 : 1/fs_orig : T_end;
xt = cos(2*pi*f0*t);

%now time to write down the Nyquist values and use the functions that we
%made

fs_below = 0.5 * fN; %bellow Nyquist
[t_sample_below, x_sample_below] = sample(t, xt, fs_below);
xrcon_below = reconstruct(t, x_sample_below, fs_below);

fs_above = 2 * fN; %above Nyquist
[t_sample_above, x_sample_above] = sample(t, xt, fs_above);
xrcon_above = reconstruct(t, x_sample_above, fs_above);

%now that we have all of our values, all that needs to be done is to plot
%everything we calculated it and copare it to the "original" cos function

figure('Name', 'Sampling and Reconstruction Demonstration');

% Plot for Below Nyquist Sampling
subplot(2,1,1);
plot(t, xt, 'b', 'LineWidth', 1);
hold on;
stem(t_sample_below, x_sample_below, 'r', 'filled', 'MarkerSize', 5);
plot(t, xrcon_below, 'g', 'LineWidth', 2);
title(sprintf('Below Nyquist (fs = %g Hz, fN = %g Hz)', fs_below, fN));
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Samples', 'Reconstructed Signal', 'Location', 'best');
grid on;

% Plot for Above Nyquist Sampling
subplot(2,1,2);
plot(t, xt, 'b', 'LineWidth', 2);
hold on;
stem(t_sample_above, x_sample_above, 'r', 'filled', 'MarkerSize', 5);
plot(t, xrcon_above, 'g', 'LineWidth', 2);
title(sprintf('Above Nyquist (fs = %g Hz, fN = %g Hz)', fs_above, fN));
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Samples', 'Reconstructed Signal', 'Location', 'best');
grid on;