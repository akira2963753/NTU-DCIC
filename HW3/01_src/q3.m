%% HW3 Q3: OFDM parameters and channel response on the subcarriers
clear; close all; clc;

%% Parameters
ch = hw3_channel();
h = ch.h;
tau = ch.tau;
Fs = ch.Fs;
L = ch.L;

p = hw_palette();

%% (a) Guard interval
Ng_min = L - 1;
[Ng, N] = ofdm_size(L);
fprintf('(a) channel taps n_i = [%s], L = %d\n', sprintf(' %d', ch.ni), L);
fprintf('    minimum guard interval = L - 1 = %d samples = %g ns\n', Ng_min, Ng_min/Fs*1e9);
fprintf('    power-of-2 candidates: 32 (< %d, too short), 64 (>= %d) -> Ng = %d = %g ns\n', ...
    Ng_min, Ng_min, Ng, Ng/Fs*1e9);

%% (b) FFT size
fprintf('\n(b) N      Ng/N     < 1/4 ?\n');
for Nc = 2.^(6:10)
    fprintf('    %4d   %.4f   %d\n', Nc, Ng/Nc, Ng/Nc < 1/4);
end
df = Fs/N;
fprintf('    -> N = %d, Ng/N = 1/%d, subcarrier spacing df = Fs/N = %.6f kHz\n', N, N/Ng, df/1e3);
fprintf('       useful symbol Tu = N/Fs = %g ns, CP = %g ns, total symbol = %d samples = %g ns\n', ...
    N/Fs*1e9, Ng/Fs*1e9, N + Ng, (N + Ng)/Fs*1e9);

%% (c) Channel response on the subcarriers
k = -N/2:N/2-1;
fk = k*df;
Hk = h*exp(-1j*2*pi*tau.'*fk);
Hk_dft = fftshift(fft(ch.hn, N));    % N-point DFT of h[n], reordered to k = -N/2 ~ N/2-1

[Hmax, imax] = max(abs(Hk));
[Hmin, imin] = min(abs(Hk));
fprintf('\n(c) f_k = k*df, k = %d ~ %d -> f_k = %.6f ~ %.6f MHz\n', k(1), k(end), fk(1)/1e6, fk(end)/1e6);
fprintf('    |H(f_k)|: max = %.4f at k = %d, min = %.4f at k = %d\n', Hmax, k(imax), Hmin, k(imin));
fprintf('    max|H(f_k) - DFT_N{h[n]}[k]| = %.1e\n', max(abs(Hk - Hk_dft)));
k_show = [-256 -255 -128 -90 -28 -1 0 1 128 255];
fprintf('       k     f_k (MHz)      |H(f_k)|   angle H(f_k) (rad)\n');
for kk = k_show
    fprintf('    %4d   %+12.6f   %.6f   %+.6f\n', kk, fk(kk + N/2 + 1)/1e6, abs(Hk(kk + N/2 + 1)), angle(Hk(kk + N/2 + 1)));
end
write_q3_table(k_show, fk(k_show + N/2 + 1), Hk(k_show + N/2 + 1));

f = linspace(-Fs/2, Fs/2, 5001);
H = h*exp(-1j*2*pi*tau.'*f);

fh = new_fig('Q3(c)', [900 700]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(f/df, abs(H), 'Color', p.muted, 'LineWidth', 0.9); hold on;
stem(k, abs(Hk), 'filled', 'Color', p.blue, 'MarkerSize', 2.2, 'LineWidth', 0.6);
set_k_axis(N);
ylim([0 2.4]);
ylabel('$|H[k]|$');
legend('$|H(f)|$, continuous', '$|H(f_k)|$', 'Location', 'northeast', 'Orientation', 'horizontal');
title('Magnitude');
nexttile;
plot(f/df, angle(H)/pi, 'Color', p.muted, 'LineWidth', 0.9); hold on;
stem(k, angle(Hk)/pi, 'filled', 'Color', p.orange, 'MarkerSize', 2.2, 'LineWidth', 0.6);
set_k_axis(N);
ylim([-1.1 1.1]);
yticks(-1:0.5:1);
ylabel('$\angle H[k]/\pi$');
title('Phase');
tl_labels(tl, sprintf('Channel response at the subcarriers: $f_k = k\\Delta f$, $\\Delta f = %.5g$ kHz, $N = %d$', df/1e3, N), ...
    'Subcarrier index $k$');
save_fig(fh, 'q3c');

%% Local functions
function write_q3_table(k_show, fk, Hk)
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, 'q3c.tex'), 'w');
fprintf(fid, '\\begin{tabular}{rrcc}\n\\toprule\n');
fprintf(fid, '$k$ & $f_k$ (MHz) & $|H[k]|$ & $\\angle H[k]$ (rad) \\\\\n\\midrule\n');
for i = 1:numel(k_show)
    fprintf(fid, '$%d$ & $%.6f$ & %.4f & $%.4f$ \\\\\n', k_show(i), fk(i)/1e6, abs(Hk(i)), angle(Hk(i)));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end

function set_k_axis(N)
xlim([-N/2 - 4, N/2 + 3]);
xticks(-N/2:64:N/2);
end
