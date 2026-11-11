%% HW4 Q2: coarse symbol timing by delay-and-correlate
clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'plot_util'));

%% Parameters
d = hw4_data();
x = d.x;
Ns = d.Ns;
p = hw_palette();
tol = 1e-12;             % relative tolerance for "equal to the maximum" (floating-point ties)

fprintf('OFDMTx: %d samples, N = %d, Ng = %d, %d symbols of %d samples\n', numel(x), d.N, d.Ng, d.nsym, Ns);
fprintf('preamble period check: ||x(33:64)-x(65:96)|| = %g, ||x(33:96)-x(97:160)|| = %g, ||CP - tail|| = %g\n', ...
    norm(x(33:64) - x(65:96)), norm(x(33:96) - x(97:160)), norm(x(1:32) - x(129:160)));

%% (a) Phi_DC with R = 64, L = 64 over the whole 640 samples
R = 64; L = 64;
phi_a = td_metric(x, R, L);
m_a = (1:numel(x)).';
[pk, ipk] = max(phi_a);
fprintf('\n(a) R = %d, L = %d: max Phi_DC = %.6f at m = %d (MATLAB max, first index)\n', R, L, pk, ipk);
plat = find(phi_a >= pk*(1 - tol));
fprintf('    samples within %.0e of the maximum: m = %s\n', tol, rng_str(plat));

fh = new_fig('Q2(a)', [900 460]);
sym_bounds(Ns, numel(x), p);
plot(m_a, phi_a, 'Color', p.blue); hold on;
plot(ipk, pk, 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 6);
text(240, 0.95*pk, sprintf('peak $(X, Y) = (%d,\\,%.4f)$', ipk, pk), 'Color', p.red, 'FontSize', 12, ...
    'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 1);
text(240, 0.93*pk, sprintf('plateau: equal values for $m = %d \\sim %d$', plat(1), plat(end)), ...
    'Color', p.ink2, 'FontSize', 11, 'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 1);
xlim([1 numel(x)]); ylim([0 1.15*pk]);
xticks([1 80:80:640]);
xlabel('$m$');
ylabel('$\Phi_{DC}(m)$');
title(sprintf('Delay-and-correlate metric, $R = %d$, $L = %d$', R, L));
save_fig(fh, 'q2a');

%% (b) First symbol only
m_b = 1:Ns;
[pk_b, ib] = max(phi_a(m_b));
plat_b = find(phi_a(m_b) >= pk_b*(1 - tol));
fprintf('\n(b) 1 <= m <= %d: m_hat_DC = argmax Phi_DC = %d (Phi_DC = %.6f)\n', Ns, ib, pk_b);
fprintf('    plateau (within %.0e of max): m = %s, %d samples, spread of values = %.3g\n', ...
    tol, rng_str(plat_b), numel(plat_b), max(phi_a(plat_b)) - min(phi_a(plat_b)));
fprintf('    periodic preamble = CP + 4 x N/4 = %d samples, R + L = %d -> plateau length %d - %d + 1 = %d\n', ...
    Ns, R + L, Ns, R + L, Ns - (R + L) + 1);
fprintf('    Phi_DC at m = %d (end of preamble) = %.6f, at m = %d = %.6f\n', Ns, phi_a(Ns), Ns + 1, phi_a(Ns + 1));

%% (c) Discard the preamble, R = 32, L = 128
z = d.z;
R = 32; L = 128;
[dc_c, nm_c] = td_metric(z, R, L);
nsym = numel(z)/Ns;
fprintf('\n(c) %d samples (%d symbols), R = %d, L = %d\n', numel(z), nsym, R, L);
[m_dc, v_dc] = per_symbol_max(dc_c, Ns);
[m_nm, v_nm] = per_symbol_max(nm_c, Ns);
for s = 1:nsym
    fprintf('    symbol %d (m = %3d ~ %3d): m_hat_DC = %3d (Phi_DC = %.6f), m_hat_NM = %3d (Phi_NM = %.6f), symbol end = %d\n', ...
        s, (s-1)*Ns + 1, s*Ns, m_dc(s), v_dc(s), m_nm(s), v_nm(s), s*Ns);
end
fprintf('    window energy sum|z(m-r)|^2 at the symbol ends: %s\n', sprintf(' %.6f', dc_c((1:nsym)*Ns)));
[~, ord] = sort(dc_c, 'descend');
fprintf('    largest Phi_DC values: %s\n', sprintf(' m=%d(%.4f)', [ord(1:6) dc_c(ord(1:6))].'));
[~, ord] = sort(nm_c, 'descend');
fprintf('    largest Phi_NM values: %s\n', sprintf(' m=%d(%.4f)', [ord(1:6) nm_c(ord(1:6))].'));

fh = new_fig('Q2(c)', [900 720]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot_metric(dc_c, m_dc, v_dc, Ns, p, '$\Phi_{DC}(m)$', '%.4f');
title(sprintf('Delay-and-correlate $\\Phi_{DC}(m)$, $R = %d$, $L = %d$', R, L));
nexttile;
plot_metric(nm_c, m_nm, v_nm, Ns, p, '$\Phi_{NM}(m)$', '%.4f');
title('Normalized metric $\Phi_{NM}(m)$');
tl_labels(tl, '', '$m$');
save_fig(fh, 'q2c');

%% (d) SNR = 12 dB
zn = d.zn;
sigma2 = d.sigma2;
[~, nm_d] = td_metric(zn, R, L);
[m_nd, v_nd] = per_symbol_max(nm_d, Ns);
fprintf('\n(d) SNR = %d dB: signal power = %.6f, noise variance = %.6f (seed %d)\n', d.snr_db, mean(abs(z).^2), sigma2, d.seed);
fprintf('    measured SNR of this realization = %.4f dB\n', 10*log10(mean(abs(z).^2)/mean(abs(zn - z).^2)));
for s = 1:nsym
    fprintf('    symbol %d: m_hat_NM = %3d (Phi_NM = %.6f), noiseless m_hat_NM = %3d (Phi_NM = %.6f), offset = %+d\n', ...
        s, m_nd(s), v_nd(s), m_nm(s), v_nm(s), m_nd(s) - m_nm(s));
end
fprintf('    Phi_NM at the true boundaries (noisy): %s\n', sprintf(' %.6f', nm_d((1:nsym)*Ns)));
ps = mean(abs(z).^2);
fprintf('    expected peak ~ (S/(S+N))^2 = %.4f\n', (ps/(ps + sigma2))^2);

fh = new_fig('Q2(d)', [900 460]);
sym_bounds(Ns, numel(z), p);
hn = plot(1:numel(z), nm_d, 'Color', p.blue); hold on;
hc = plot(1:numel(z), nm_c, ':', 'Color', p.muted, 'LineWidth', 1.2);
xlim([1 numel(z)]); ylim([0 1.25]); xticks([1 40:40:480]);
peak_labels(m_nd, v_nd, p.red);
legend([hn hc], {sprintf('SNR = %d dB', d.snr_db), 'noiseless'}, 'Location', 'north', 'Orientation', 'horizontal');
xlabel('$m$');
ylabel('$\Phi_{NM}(m)$');
title(sprintf('Normalized metric with noise, SNR = %d dB, $R = %d$, $L = %d$', d.snr_db, R, L));
save_fig(fh, 'q2d');

%% Local functions
function [mi, v] = per_symbol_max(phi, Ns)
% argmax of phi within every Ns-sample block (index relative to the whole record)
nsym = numel(phi)/Ns;
mi = zeros(nsym, 1);
v = zeros(nsym, 1);
for s = 1:nsym
    [v(s), k] = max(phi((s-1)*Ns + (1:Ns)));
    mi(s) = (s-1)*Ns + k;
end
end

function plot_metric(phi, mi, v, Ns, p, yl, fmt)
sym_bounds(Ns, numel(phi), p);
plot(1:numel(phi), phi, 'Color', p.blue); hold on;
xlim([1 numel(phi)]); ylim([0 1.2*max(v)]);
xticks([1 40:40:numel(phi)]);
peak_labels(mi, v, p.red, fmt);
ylabel(yl);
end

function sym_bounds(Ns, M, p)
for b = Ns:Ns:M
    xline(b, '--', 'Color', p.muted, 'LineWidth', 0.8, 'HandleVisibility', 'off', 'Layer', 'bottom'); hold on;
end
end

function s = rng_str(k)
% compact text for a list of indices
if isempty(k), s = '[]'; return; end
if all(diff(k) == 1), s = sprintf('%d ~ %d', k(1), k(end)); else, s = sprintf(' %d', k); end
end
