% MATLAB Script for OOK Modulation Simulation and Plotting (Optimized without tightlayout)

% Create 'Pictures' folder if it doesn't exist
if ~exist('Pictures', 'dir')
    mkdir('Pictures');
    disp('Created "Pictures" folder.');
end

% --- Parameters ---
Rb = 1000; % Bit rate (bits/second)
T = 1/Rb; % Bit duration
fc = 20e3; % Carrier frequency (e.g., 20 kHz)
Ac = 1; % Carrier amplitude

fs_sampling_factor = 200; % Factor for sampling rate (fs_sampling = fs_sampling_factor * Rb)
fs_sampling = fs_sampling_factor * Rb; % High sampling frequency for smooth plots
Ts_sampling = 1/fs_sampling;

num_bits = 10; % Number of bits to simulate
bits = randi([0, 1], 1, num_bits); % Random binary data (0s and 1s)

t = 0:Ts_sampling:(num_bits * T - Ts_sampling);

% --- OOK Baseband Signal (m(t)) ---
m_t_ook = zeros(size(t));
for i = 1:num_bits
    idx_start = round((i-1)*T/Ts_sampling) + 1;
    idx_end = round(i*T/Ts_sampling);
    if idx_start > length(t) || idx_end > length(t)
        idx_end = length(t); 
    end
    if idx_start > idx_end 
        continue; 
    end
    m_t_ook(idx_start : idx_end) = bits(i);
end

% --- OOK Complex Envelope (g(t)) ---
g_t_ook = Ac * m_t_ook; % g(t) is real and positive for OOK

% --- Fourier Transform of g(t) (G(f)) ---
N_fft_g = length(g_t_ook); 
N_fft_g = 2^nextpow2(N_fft_g); % Pad to nearest power of 2
g_t_ook_padded_for_fft = [g_t_ook, zeros(1, N_fft_g - length(g_t_ook))];

G_f_ook = fftshift(fft(g_t_ook_padded_for_fft)) / fs_sampling;
f_g_ook = (-N_fft_g/2 : N_fft_g/2 - 1) * (fs_sampling / N_fft_g);

% --- OOK Passband Signal (s(t)) ---
s_t_ook = Ac * m_t_ook .* cos(2*pi*fc*t); % OOK Passband Signal: s(t) = A_c * m(t) * cos(2*pi*f_c*t)

% --- Fourier Transform of s(t) (S(f)) ---
N_fft_s = length(s_t_ook);
N_fft_s = 2^nextpow2(N_fft_s); 
s_t_ook_padded_for_fft = [s_t_ook, zeros(1, N_fft_s - length(s_t_ook))];

S_f_ook = fftshift(fft(s_t_ook_padded_for_fft)) / fs_sampling;
f_s_ook = (-N_fft_s/2 : N_fft_s/2 - 1) * (fs_sampling / N_fft_s);

% --- Plotting and Exporting ---

% Plot 1: m(t) and g(t) in time domain
fig1 = figure('Units', 'inches', 'Position', [0 0 8 6]); % Consistent figure size

% Get subplot handles for manual positioning
ax1 = subplot(2,1,1);
plot(t, m_t_ook, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); % Blue
title('Señal Modulante m(t) (Tren de Pulsos NRZ Unipolar)', 'FontSize', 11);
xlabel('Tiempo (s)', 'FontSize', 9);
ylabel('Amplitud', 'FontSize', 9);
ylim([-0.2 1.2]);
grid on;
set(ax1, 'FontSize', 9); % Apply font size to axes

% Adjust position of the first subplot to make space for sgtitle
ax1.Position = ax1.Position + [0 0.05 0 0]; % Move up slightly

ax2 = subplot(2,1,2);
plot(t, real(g_t_ook), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % Orange
title('Envolvente Compleja g(t) para OOK', 'FontSize', 11);
xlabel('Tiempo (s)', 'FontSize', 9);
ylabel('Amplitud', 'FontSize', 9);
ylim([-0.2 1.2]);
grid on;
set(ax2, 'FontSize', 9); % Apply font size to axes

% Adjust position of the second subplot to move down slightly
ax2.Position = ax2.Position - [0 0.05 0 0]; 

sgtitle('Envolvente Compleja OOK en el Dominio del Tiempo', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig1, 'Pictures/ook_mt_gt.png', 'Resolution', 400); 
disp('Saved Pictures/ook_mt_gt.png');

% Plot 2: G(f) spectrum (baseband)
fig2 = figure('Units', 'inches', 'Position', [0 0 8 4]); % Consistent figure size
plot(f_g_ook, abs(G_f_ook), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]); % Purple
title('Espectro de la Envolvente Compleja G(f) para OOK', 'FontSize', 11);
xlabel('Frecuencia (Hz)', 'FontSize', 9);
ylabel('|G(f)|', 'FontSize', 9, 'Interpreter', 'latex');
xlim([-3*Rb 3*Rb]); 
grid on;
set(gca, 'FontSize', 9);

sgtitle('Espectro de la Envolvente Compleja OOK', 'FontSize', 13, 'FontWeight', 'bold');
% No tightlayout needed for single subplot figures, but adjust Position for more margin
ax = gca; % Get current axes handle
ax.Position = ax.Position + [0.02 0.05 -0.04 -0.08]; % [left bottom width height] -> adjust to make more room
exportgraphics(fig2, 'Pictures/ook_gf.png', 'Resolution', 400); 
disp('Saved Pictures/ook_gf.png');

% Plot 3: S(f) spectrum (passband)
fig3 = figure('Units', 'inches', 'Position', [0 0 8 4]); % Consistent figure size
plot(f_s_ook / 1e3, abs(S_f_ook), 'LineWidth', 2, 'Color', [0 0.6 0.2]); % Dark Green
title('Espectro de Señal Modulada OOK (Pasabanda)', 'FontSize', 11);
xlabel('Frecuencia (kHz)', 'FontSize', 9);
ylabel('Magnitud Espectral', 'FontSize', 9);
xlim([ (fc - 2*Rb)/1e3, (fc + 2*Rb)/1e3]); 
grid on;
set(gca, 'FontSize', 9);

sgtitle('Espectro de la Señal OOK Modulada', 'FontSize', 13, 'FontWeight', 'bold');
% No tightlayout needed for single subplot figures, but adjust Position for more margin
ax = gca;
ax.Position = ax.Position + [0.02 0.05 -0.04 -0.08]; % [left bottom width height] -> adjust to make more room
exportgraphics(fig3, 'Pictures/ook_passband_spectrum.png', 'Resolution', 400); 
disp('Saved Pictures/ook_passband_spectrum.png');

close all;
disp('OOK plots generated and saved in the "Pictures" folder.');
