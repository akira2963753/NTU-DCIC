%% HW3 Q4: OFDM symbol through the four-path channel
clear; close all; clc;

%% Parameters
ch = hw3_channel();    % seeds rng(39) and draws gamma_i first; the QPSK bits below continue that stream
h = ch.h;
hn = ch.hn;
tau = ch.tau;
Fs = ch.Fs;
L = ch.L;
[Ng, N] = ofdm_size(L);

p = hw_palette();
k = -N/2:N/2-1;
xl_k = 'Subcarrier index $k$';
xl_n = 'Sample index $n$';
fprintf('Fs = %g MHz, L = %d, Ng = %d, N = %d\n', Fs/1e6, L, Ng, N);

%% (a) QPSK data on all N subcarriers, IFFT and cyclic prefix
% Black points of HW1 Fig. 1: (0,1)->45, (1,1)->135, (1,0)->225, (0,0)->315 deg
b = randi([0 1], 2, N);
X = ((1 - 2*b(1, :)) + 1j*(2*b(2, :) - 1))/sqrt(2);    % X(k + N/2 + 1) = X[k]

% Unitary DFT pair: x[n] = 1/sqrt(N) sum_k X[k] e^{j2pi kn/N}
x = sqrt(N)*ifft(ifftshift(X));
xcp = [x(N-Ng+1:N), x];
ncp = 0:N+Ng-1;

fprintf('\n(a) first subcarriers:\n');
for kk = k(1:4)
    fprintf('    k = %4d: (b0,b1) = (%d,%d), X[k] = %+.4f%+.4fj\n', kk, b(:, kk + N/2 + 1), ...
        real(X(kk + N/2 + 1)), imag(X(kk + N/2 + 1)));
end
fprintf('    symbol counts (45/135/225/315 deg) = %s\n', ...
    sprintf('%d ', histcounts(mod(round(rad2deg(angle(X))), 360), [0 90 180 270 360])));
fprintf('    x[n]: length %d + CP %d = %d samples, mean|x|^2 = %.6f, max|x| = %.4f\n', ...
    N, Ng, numel(xcp), mean(abs(x).^2), max(abs(x)));
fprintf('    CP check: max|x_cp[0..%d] - x[%d..%d]| = %.1e\n', Ng - 1, N - Ng, N - 1, ...
    max(abs(xcp(1:Ng) - xcp(N+1:N+Ng))));

fh = new_fig('Q4(a) X[k]', [900 900]);
tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1 2]);
stem_qpsk(k, real(X), p.blue);
set_k_axis(N);
ylabel('$\mathrm{Re}\{X[k]\}$');
title('Real part, all $N = 512$ subcarriers');
nexttile([1 2]);
stem_qpsk(k, imag(X), p.orange);
set_k_axis(N);
ylabel('$\mathrm{Im}\{X[k]\}$');
title('Imaginary part, all $N = 512$ subcarriers');
kz = -16:15;    % zoom so the individual symbols are readable
nexttile;
stem_qpsk(kz, real(X(kz + N/2 + 1)), p.blue);
xlim([kz(1) - 1, kz(end) + 1]); xticks(-16:8:16);
ylabel('$\mathrm{Re}\{X[k]\}$');
title('Real part, zoom $k = -16 \sim 15$');
nexttile;
stem_qpsk(kz, imag(X(kz + N/2 + 1)), p.orange);
xlim([kz(1) - 1, kz(end) + 1]); xticks(-16:8:16);
ylabel('$\mathrm{Im}\{X[k]\}$');
title('Imaginary part, zoom $k = -16 \sim 15$');
tl_labels(tl, 'Frequency-domain QPSK data $X[k]$ (black points of HW1)', xl_k);
save_fig(fh, 'q4a1');

fh = new_fig('Q4(a) x[n]', [900 640]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot_cp_signal(ncp, real(xcp), Ng, N, p, p.blue);
ylabel('$\mathrm{Re}\{x[n]\}$');
legend(sprintf('CP, $0 \\le n \\le %d$', Ng - 1), sprintf('Copied from $%d \\le n \\le %d$', N, N + Ng - 1), ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
nexttile;
plot_cp_signal(ncp, imag(xcp), Ng, N, p, p.orange);
ylabel('$\mathrm{Im}\{x[n]\}$');
tl_labels(tl, sprintf('Time-domain OFDM symbol with cyclic prefix, $N + N_g = %d$ samples', N + Ng), xl_n);
save_fig(fh, 'q4a2');

%% (b) Linear convolution with h[n], CP removal and DFT
y = conv(hn, xcp);
ny = 0:numel(y)-1;
yd = y(Ng+1:Ng+N);
Y = fftshift(fft(yd))/sqrt(N);

fprintf('\n(b) y[n] = h[n] (*) x[n]: length %d + %d - 1 = %d samples\n', L, N + Ng, numel(y));
fprintf('    kept samples n = %d ~ %d after removing the CP\n', Ng, Ng + N - 1);
fprintf('    max|y_kept - cconv(h, x, N)| = %.1e\n', max(abs(yd - cconv(hn, x, N))));
fprintf('    |Y[k]|: min = %.4f, max = %.4f, mean|Y|^2 = %.6f\n', min(abs(Y)), max(abs(Y)), mean(abs(Y).^2));

fh = new_fig('Q4(b) time', [900 640]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot_rx_signal(ny, real(y), Ng, N, L, p, p.blue);
ylabel('$\mathrm{Re}\{y[n]\}$');
legend('CP (removed)', 'Channel tail (removed)', 'Location', 'northoutside', 'Orientation', 'horizontal');
nexttile;
plot_rx_signal(ny, imag(y), Ng, N, L, p, p.orange);
ylabel('$\mathrm{Im}\{y[n]\}$');
tl_labels(tl, sprintf('Received signal $y[n] = h[n] \\otimes x[n]$, $%d$ samples', numel(y)), xl_n);
save_fig(fh, 'q4b1');

fh = new_fig('Q4(b) freq', [900 420]);
stem(k, abs(Y), 'filled', 'Color', p.blue, 'MarkerSize', 2.2, 'LineWidth', 0.6); hold on;
yline(1, '--', 'Color', p.red, 'LineWidth', 1.1);
set_k_axis(N);
ylim([0 2.4]);
xlabel(xl_k);
ylabel('$|Y[k]|$');
legend('$|Y[k]|$', '$|X[k]| = 1$', 'Location', 'northeast', 'Orientation', 'horizontal');
title('Received frequency-domain signal after CP removal and $N$-point DFT');
save_fig(fh, 'q4b2');

%% (c) Channel estimate H[k] = Y[k]/X[k] vs. Q3(c)
Hest = Y./X;
Hq3 = h*exp(-1j*2*pi*tau.'*(k*Fs/N));    % H(f_k) of Q3(c)
err = abs(Hest - Hq3);
fprintf('\n(c) max|H_est[k] - H(f_k)| = %.1e, max||H_est| - |H(f_k)|| = %.1e\n', ...
    max(err), max(abs(abs(Hest) - abs(Hq3))));

% Same chain with other CP lengths, to show where the match breaks
Ng_try = [32 49 50 64];
err_try = zeros(numel(Ng_try), N);
for i = 1:numel(Ng_try)
    err_try(i, :) = abs(ofdm_chain(X, hn, Ng_try(i))./X - Hq3);
    fprintf('    Ng = %2d: max|H_est - H(f_k)| = %.2e, mean = %.2e\n', Ng_try(i), max(err_try(i, :)), mean(err_try(i, :)));
end

fh = new_fig('Q4(c)', [900 980]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem(k, abs(Hest), 'filled', 'Color', p.blue, 'MarkerSize', 2.2, 'LineWidth', 0.6); hold on;
plot(k(1:4:end), abs(Hq3(1:4:end)), 'o', 'Color', p.ink, 'MarkerSize', 4, 'LineWidth', 0.7);
set_k_axis(N);
ylim([0 2.4]);
ylabel('$|H[k]|$');
legend('$|Y[k]/X[k]|$', '$|H(f_k)|$ of Q3(c), every 4th $k$', 'Location', 'northeast', 'Orientation', 'horizontal');
title('Magnitude');
nexttile;
stem(k, angle(Hest)/pi, 'filled', 'Color', p.orange, 'MarkerSize', 2.2, 'LineWidth', 0.6); hold on;
plot(k(1:4:end), angle(Hq3(1:4:end))/pi, 'o', 'Color', p.ink, 'MarkerSize', 4, 'LineWidth', 0.7);
set_k_axis(N);
ylim([-1.1 1.1]);
yticks(-1:0.5:1);
ylabel('$\angle H[k]/\pi$');
title('Phase');
nexttile;
% Ng = 50 and 64 give identical results, so they share one marker set
cols = {p.red, p.muted, p.blue};
for i = 1:3
    semilogy(k, max(err_try(i, :), 1e-17), '.', 'Color', cols{i}, 'MarkerSize', 6); hold on;
end
set_k_axis(N);
ylim([1e-17 10]);
yticks(10.^(-16:4:0));
ylabel('$|H[k] - H(f_k)|$');
legend('$N_g = 32$', '$N_g = 49$', '$N_g = 50$ and $64$ (identical)', 'Location', 'east');
title('Estimation error vs.\ CP length (channel memory $L - 1 = 50$)');
tl_labels(tl, sprintf('Estimated channel $H[k] = Y[k]/X[k]$ ($N_g = %d$) vs.\\ Q3(c)', Ng), xl_k);
save_fig(fh, 'q4c');

%% Local functions
function Y = ofdm_chain(X, hn, Ng)
% IFFT -> add CP -> linear convolution -> remove CP -> FFT
N = numel(X);
x = sqrt(N)*ifft(ifftshift(X));
y = conv(hn, [x(N-Ng+1:N), x]);
Y = fftshift(fft(y(Ng+1:Ng+N)))/sqrt(N);
end

function stem_qpsk(k, v, color)
stem(k, v, 'filled', 'Color', color, 'MarkerSize', 2.6, 'LineWidth', 0.7);
ylim([-1 1]);
yticks([-1/sqrt(2) 0 1/sqrt(2)]);
yticklabels({'$-1/\sqrt{2}$', '$0$', '$1/\sqrt{2}$'});
end

function plot_cp_signal(n, v, Ng, N, p, color)
xregion(-0.5, Ng - 0.5, 'FaceColor', p.muted, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
xregion(N - 0.5, N + Ng - 0.5, 'FaceColor', p.aqua, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hold on;
plot(n, v, 'Color', color, 'LineWidth', 0.9, 'HandleVisibility', 'off');
xlim([0 N + Ng]);
xticks([0 Ng 128:128:N N + Ng]);
ylim([-3 3]);
end

function plot_rx_signal(n, v, Ng, N, L, p, color)
xregion(-0.5, Ng - 0.5, 'FaceColor', p.muted, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
xregion(N + Ng - 0.5, N + Ng + L - 1.5, 'FaceColor', p.red, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
hold on;
plot(n, v, 'Color', color, 'LineWidth', 0.9, 'HandleVisibility', 'off');
xlim([0 numel(n)]);
xticks([0 Ng 128:128:N N + Ng N + Ng + L - 1]);
ylim([-3 3]);
end

function set_k_axis(N)
xlim([-N/2 - 4, N/2 + 3]);
xticks(-N/2:64:N/2);
end
