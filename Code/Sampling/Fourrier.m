function [xhat, ck] = ffs(xt, t, n, T)

    %first we need to limit the t values to -T/2=< t <= T/2
    t_period_indices = find(t >= -T/2 & t <= T/2);
    t_period = t(t_period_indices);
    xt_period = xt(t_period_indices);
    
    %making the empty vector
    ck = zeros(1, 2*n + 1);
    
    %we need to loop from -n to n since that's the limits of the sum
    k_values = -n:n;
    for idx = 1:(2*n + 1)
        k = k_values(idx);
        y = xt_period .* exp(-1j * 2 * pi * k / T * t_period);
        %here we will use the trapezoidal rule for numerical integration
        %which is a function in matlab where t_period will represent
        %the x coordinates and y represents the y coordinates
        integral_val = trapz(t_period, y);
        ck(idx) = (1/T) * integral_val;
    end
    
    %now to compute the sum, first w eneed an empty list
    xhat = zeros(size(t));
    for idx = 1:(2*n + 1)
        k = k_values(idx);
        %now we simply add the previous value of xhat and multiply ck which we
        %calculated previously with the rest of the function
        xhat = xhat + ck(idx) * exp(1j * 2 * pi * k / T * t);
    end
end

%Now to define our signal
fs = 1000;
t_max = 10;
t = -t_max : 1/fs : t_max; %Time range from -t_max to +t_max

%now to make our sinc function
A = 1;
B = 1;
%I added these two variables in order to be able to easily manipulate the
%amplitude and the bandwith factor
x_t = A * sinc(B*t);

%part a and b)

T_fixed = 15;
n_vals = [1, 5, 15, 50];

figure('Name', 'Approximation with Varying n');
plot(t, x_t, 'g', 'LineWidth', 2.5, 'DisplayName', 'Original x(t)');
hold on;
grid on;
colors = lines(length(n_vals)); %Used this function in order to automatically
%allocate colors for each graph which was followed by a loop to generate
%each graph, making the code more compact and better looking

for i = 1:length(n_vals)
    n_current = n_vals(i);
    [xhat, ck] = ffs(x_t, t, n_current, T_fixed);
    xhat = real(xhat);
    plot(t, xhat, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', ['n = ' num2str(n_current)]);
end

title(['Fourier Series Approximation for Fixed T = ' num2str(T_fixed)]);
xlabel('Time (t)');
ylabel('Amplitude');
legend('show', 'Location', 'northeast');
xlim([-t_max, t_max]); %for the graphs I'll be limiting it
%to -1-T/2 and 1+T/2 since any additional information is not needed
hold off;

%part c)

n_fixed = 30; %large number of n for best accuracy
T_vals = [6, 15, 30, 80];

figure('Name', 'Approximation with Varying T');
plot(t, x_t, 'g', 'LineWidth', 2.5, 'DisplayName', 'Original x(t)');
hold on;
grid on;
colors = lines(length(T_vals));

for i = 1:length(T_vals)
    T_current = T_vals(i);
    [xhat, ck] = ffs(x_t, t, n_fixed, T_current);
    xhat = real(xhat);
    plot(t, xhat, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', ['T = ' num2str(T_current)]);
end

title(['Fourier Series Approximation for Fixed n = ' num2str(n_fixed)]);
xlabel('Time (t)');
ylabel('Amplitude');
legend('show', 'Location', 'northeast');
xlim([-t_max, t_max]);
hold off;

%part d)

n_range = 1:50; % Range of n values to analyze
error_vs_n = zeros(size(n_range));
T_fixed_error = 20;

for i = 1:length(n_range)
    n_current = n_range(i);
    [xhat, ck] = ffs(x_t, t, n_current, T_fixed_error);
    % Calculate error energy: E = integral(|x(t) - xhat(t)|^2) dt
    error_y = abs(x_t - xhat).^2;
    %again, using the trapezoid rule to calculate integral where t is the x
    %representation and error_y is the y
    error_vs_n(i) = trapz(t, error_y);
end

figure('Name', 'Error vs. n');
semilogy(n_range, error_vs_n, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
grid on;
title(['Approximation Error Energy vs. n (for T = ' num2str(T_fixed_error) ')']);
xlabel('Number of Harmonics (n)');
ylabel('Error Energy E(n, T)  (log scale)');

%part e)

n_fixed_error = 30; %50 is sufficiently large
T_range = 1:0.5:15;  % Range of T values to analyze
error_vs_T = zeros(size(T_range));

for i = 1:length(T_range)
    T_current = T_range(i);
    [xhat, ck] = ffs(x_t, t, n_fixed_error, T_current);
    error_integrand = abs(x_t - xhat).^2;
    error_vs_T(i) = trapz(t, error_integrand);
end

figure('Name', 'Error vs. T');
plot(T_range, error_vs_T, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3);
grid on;
title(['Approximation Error Energy vs. T (for n = ' num2str(n_fixed_error) ')']);
xlabel('Period (T)');
ylabel('Error Energy E(n, T)');