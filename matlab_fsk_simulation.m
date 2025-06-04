% MATLAB Script for FSK Modulation Simulation and Plotting (Optimized without tightlayout)

% Create 'Pictures' folder if it doesn't exist
if ~exist('Pictures', 'dir')
    mkdir('Pictures');
    disp('Created "Pictures" folder.');
end

% --- Parameters ---
Rb = 1000; % Bit rate (bits/second)
T = 1/Rb; % Bit duration

fc_common = 20e3; % Common carrier frequency (e.g., 20 kHz)
f_dev_fsk = 5e3; % Frequency deviation for FSK (e.g., 5 kHz from center)
fc_mark = fc_common + f_dev_fsk; % Frequency for '1' (e.g., 20k + 5k = 25 kHz)
fc_space = fc_common - f_dev_fsk; % Frequency for '0' (e.g., 20k - 5k = 15 kHz)
Ac = 1; % Carrier amplitude

fs_sampling_factor = 200; % Factor for sampling rate
fs_sampling = fs_sampling_factor * max(fc_mark, fc_space); % Ensure high enough sampling
Ts_sampling = 1/fs_sampling;

num_bits = 4; % Reduced number of bits for clearer time-domain visualization of FSK changes
bits = randi([0, 1], 1, num_bits); % Random binary data (0s and 1s)

t = 0:Ts_sampling:(num_bits * T - Ts_sampling);
t_samples_per_bit = floor(T/Ts_sampling);

% --- FSK Complex Envelope (g(t)) ---
g_t_fsk = zeros(size(t));
for i = 1:num_bits
    idx_start = round((i-1)*T/Ts_sampling) + 1;
    idx_end = round(i*T/Ts_sampling);
    if idx_start > length(t) || idx_end > length(t)
        idx_end = length(t);
    end
    if idx_start > idx_end
        continue;
    end
    t_segment = t(idx_start : idx_end);
    
    if bits(i) == 1
        f_relative = (fc_mark - fc_common); 
    else
        f_relative = (fc_space - fc_common); 
    end
    g_t_fsk(idx_start : idx_end) = Ac * exp(1j * 2 * pi * f_relative * t_segment);
end
g_t_fsk = g_t_fsk(1:length(t)); 

% --- Fourier Transform of g(t) (G(f)) ---
N_fft_g = length(g_t_fsk);
N_fft_g = 2^nextpow2(N_fft_g); 
g_t_fsk_padded_for_fft = [g_t_fsk, zeros(1, N_fft_g - length(g_t_fsk))];

G_f_fsk = fftshift(fft(g_t_fsk_padded_for_fft)) / fs_sampling;
f_g_fsk = (-N_fft_g/2 : N_fft_g/2 - 1) * (fs_sampling / N_fft_g);

% --- FSK Passband Signal (s(t)) ---
s_t_fsk = zeros(size(t));
for i = 1:num_bits
    idx_start = round((i-1)*T/Ts_sampling) + 1;
    idx_end = round(i*T/Ts_sampling);
    if idx_start > length(t) || idx_end > length(t)
        idx_end = length(t); 
    end
    if idx_start > idx_end 
        continue; 
    end
    t_segment = t(idx_start : idx_end);
    
    if bits(i) == 1
        s_t_fsk(idx_start : idx_end) = Ac * cos(2*pi*fc_mark*t_segment);
    else
        s_t_fsk(idx_start : idx_end) = Ac * cos(2*pi*fc_space*t_segment);
    end
end

% --- Fourier Transform of s(t) (S(f)) ---
N_fft_s = length(s_t_fsk);
N_fft_s = 2^nextpow2(N_fft_s); 
s_t_fsk_padded_for_fft = [s_t_fsk, zeros(1, N_fft_s - length(s_t_fsk))];

S_f_fsk = fftshift(fft(s_t_fsk_padded_for_fft)) / fs_sampling;
f_s_fsk = (-N_fft_s/2 : N_fft_s/2 - 1) * (fs_sampling / N_fft_s);

% --- Plotting and Exporting ---

% Plot 1: Real and Imaginary parts of g(t)
fig1 = figure('Units', 'inches', 'Position', [0 0 8 6]);
ax1 = subplot(2,1,1);
plot(t, real(g_t_fsk), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); % Blue
hold on;
plot(t, imag(g_t_fsk), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % Orange
hold off;
title('Componentes de la Envolvente Compleja g(t) para FSK', 'FontSize', 11);
xlabel('Tiempo (s)', 'FontSize', 9);
ylabel('Amplitud', 'FontSize', 9);
grid on;
legend('Re\{g(t)\}', 'Im\{g(t)\}', 'FontSize', 8, 'Location', 'best', 'Interpreter', 'latex');
set(ax1, 'FontSize', 9); % Apply font size to axes
ax1.Position = ax1.Position + [0 0.05 0 0]; % Move up slightly

ax2 = subplot(2,1,2);
plot(t, abs(g_t_fsk), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]); % Purple
title('Magnitud de la Envolvente Compleja |g(t)| para FSK', 'FontSize', 11);
xlabel('Tiempo (s)', 'FontSize', 9);
ylabel('Magnitud', 'FontSize', 9);
ylim([-0.2 1.2]); % Should be constant (Ac=1)
grid on;
set(ax2, 'FontSize', 9); % Apply font size to axes
ax2.Position = ax2.Position - [0 0.05 0 0]; % Move down slightly

sgtitle('Envolvente Compleja FSK en el Dominio del Tiempo', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig1, 'Pictures/fsk_gt_components.png', 'Resolution', 400); 
disp('Saved Pictures/fsk_gt_components.png');

% Plot 3: Magnitude of g(t) (separate for report structure, if needed)
fig2 = figure('Units', 'inches', 'Position', [0 0 8 4]);
ax = gca;
plot(t, abs(g_t_fsk), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]); % Purple
title('Magnitud de la Envolvente Compleja |g(t)| para FSK', 'FontSize', 11);
xlabel('Tiempo (s)', 'FontSize', 9);
ylabel('Magnitud', 'FontSize', 9);
ylim([-0.2 1.2]);
grid on;
set(ax, 'FontSize', 9);
ax.Position = ax.Position + [0.02 0.05 -0.04 -0.08]; % Adjust position for more margin
sgtitle('Magnitud de la Envolvente Compleja FSK', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig2, 'Pictures/fsk_gt_magnitude.png', 'Resolution', 400); 
disp('Saved Pictures/fsk_gt_magnitude.png');


% Plot 4: G(f) spectrum (baseband)
fig3 = figure('Units', 'inches', 'Position', [0 0 8 4]);
ax = gca;
plot(f_g_fsk, abs(G_f_fsk), 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880]); % Green
title('Espectro de la Envolvente Compleja G(f) para FSK', 'FontSize', 11);
xlabel('Frecuencia (Hz)', 'FontSize', 9);
ylabel('|G(f)|', 'FontSize', 9, 'Interpreter', 'latex');
xline(-(fc_mark - fc_common), 'r--', '$\mathbf{-\Delta f}$', 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left', 'FontSize', 9, 'Interpreter', 'latex');
xline( (fc_mark - fc_common), 'r--', '$\mathbf{+\Delta f}$', 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right', 'FontSize', 9, 'Interpreter', 'latex');
xlim([-3*max(f_dev_fsk, Rb) 3*max(f_dev_fsk, Rb)]); 
grid on;
set(ax, 'FontSize', 9);
ax.Position = ax.Position + [0.02 0.05 -0.04 -0.08]; % Adjust position for more margin
sgtitle('Espectro de la Envolvente Compleja FSK', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig3, 'Pictures/fsk_gf.png', 'Resolution', 400); 
disp('Saved Pictures/fsk_gf.png');

% Plot 5: S(f) spectrum (passband)
fig4 = figure('Units', 'inches', 'Position', [0 0 8 4]);
ax = gca;
plot(f_s_fsk / 1e3, abs(S_f_fsk), 'LineWidth', 2, 'Color', [0 0.4 0.6]); % Dark Cyan
title('Espectro de Señal Modulada FSK (Pasabanda)', 'FontSize', 11);
xlabel('Frecuencia (kHz)', 'FontSize', 9);
ylabel('Magnitud Espectral', 'FontSize', 9);
xline(fc_space/1e3, 'r--', '$f_{\text{space}}$', 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left', 'FontSize', 9, 'Interpreter', 'latex');
xline(fc_mark/1e3, 'r--', '$f_{\text{mark}}$', 'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right', 'FontSize', 9, 'Interpreter', 'latex');
xlim([ (min(fc_mark, fc_space) - 2*Rb)/1e3, (max(fc_mark, fc_space) + 2*Rb)/1e3]); 
grid on;
set(ax, 'FontSize', 9);
ax.Position = ax.Position + [0.02 0.05 -0.04 -0.08]; % Adjust position for more margin
sgtitle('Espectro de la Señal FSK Modulada', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig4, 'Pictures/fsk_passband_spectrum.png', 'Resolution', 400); 
disp('Saved Pictures/fsk_passband_spectrum.png');

close all;
disp('FSK plots generated and saved in the "Pictures" folder.');
