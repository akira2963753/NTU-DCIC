%% HW1 Q1: single-carrier QPSK / O-QPSK
clear; close all; clc;

%% Parameters
seed = 39;
Tsym = 0.1e-6;
fc = 40e6;
fs = 640e6;
sps = 64;
nsym = 16;
pb_gain = 1;    % passband: x_PB(t) = Re{x_BB(t) e^(j2*pi*fc*t)}, course slides Ch.2 p.16

p = hw_palette();
N = sps*nsym;
t = (0:N-1)/fs;
t_us = t*1e6;
Tsym_us = Tsym*1e6;
xl_t = '$t$ ($\mu$s)';

%% Bits and QPSK mapping (black points in Fig. 1)
b = gen_bits(seed);
fprintf('b_0 ~ b_31 = %s\n', sprintf('%d', b));

% (0,1)->45, (1,1)->135, (1,0)->225, (0,0)->315 deg
b0 = b(1:2:end);
b1 = b(2:2:end);
a = ((1 - 2*b0) + 1j*(2*b1 - 1))/sqrt(2);

fprintf('\n k  (b2k,b2k+1)  a_k\n');
for k = 0:nsym-1
    fprintf('%2d     (%d,%d)     %+.4f%+.4fj  (%4.0f deg)\n', k, b0(k+1), b1(k+1), ...
        real(a(k+1)), imag(a(k+1)), rad2deg(angle(a(k+1))));
end

%% (a) Baseband QPSK
xbb = repelem(a, sps);

fh = new_fig('Q1(a)');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_us, real(xbb), 'Color', p.blue);
sym_lines(0:nsym, Tsym_us, p.muted);
set_amp_axis([0 nsym*Tsym_us]);
ylabel('$\mathrm{Re}\{x_{BB}(t)\}$');
nexttile;
plot(t_us, imag(xbb), 'Color', p.orange);
sym_lines(0:nsym, Tsym_us, p.muted);
set_amp_axis([0 nsym*Tsym_us]);
ylabel('$\mathrm{Im}\{x_{BB}(t)\}$');
tl_labels(tl, 'QPSK baseband signal', xl_t);
save_fig(fh, 'q1a');

%% (b) Passband signal, first 8 symbols
xpb = pb_gain*real(xbb.*exp(1j*2*pi*fc*t));
n8 = 8*sps;

fh = new_fig('Q1(b)', [900 420]);
plot(t_us(1:n8), xpb(1:n8), 'Color', p.aqua, 'LineWidth', 1.3);
xline((0:8)*Tsym_us, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');    % spec: black dashed boundaries
set_amp_axis([0 8*Tsym_us]);
grid off; set(gca, 'YGrid', 'on');    % x grid coincides with the boundaries
xlabel(xl_t); ylabel('$x_{PB}(t)$');
title('QPSK passband signal (first 8 symbols)');
save_fig(fh, 'q1b');

%% (c) Baseband vs passband, first 2 symbols
n2 = 2*sps;

fh = new_fig('Q1(c)', [900 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile;
plot(t_us(1:n2), real(xbb(1:n2)), 'Color', p.blue); hold on;
plot(t_us(1:n2), imag(xbb(1:n2)), '--', 'Color', p.orange);
sym_lines(0:2, Tsym_us, p.muted);
set_amp_axis([0 2*Tsym_us]);
ylabel('Amplitude');
legend('$\mathrm{Re}\{x_{BB}\}$', '$\mathrm{Im}\{x_{BB}\}$', 'Location', 'eastoutside');
title('Baseband');
ax2 = nexttile;
plot(t_us(1:n2), xpb(1:n2), '.-', 'Color', p.aqua, 'MarkerSize', 9, 'LineWidth', 1.2);
sym_lines(0:2, Tsym_us, p.muted);
set_amp_axis([0 2*Tsym_us]);
ylabel('Amplitude');
legend('$x_{PB}$', 'Location', 'eastoutside');
title(sprintf('Passband ($f_c$ = %g MHz, %d samples per carrier cycle)', fc/1e6, fs/fc));
linkaxes([ax1 ax2], 'x');
tl_labels(tl, '', xl_t);
save_fig(fh, 'q1c');

%% (d) Magnitude spectra
nfft = 8*N;    % zero padding for a smoother spectrum
f_MHz = (-nfft/2:nfft/2-1)*fs/nfft/1e6;
Xbb = fftshift(fft(xbb, nfft));
Xpb = fftshift(fft(xpb, nfft));

fh = new_fig('Q1(d)', [900 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(f_MHz, mag_db(Xbb), 'Color', p.blue, 'LineWidth', 1.1);
set_db_axis();
ylabel('$|X_{BB}(f)|$ (dB)');
title('Baseband spectrum');
nexttile;
plot(f_MHz, mag_db(Xpb), 'Color', p.aqua, 'LineWidth', 1.1);
fc_line(-fc/1e6, '$-f_c$', 'left', p);
fc_line(fc/1e6, '$f_c$', 'right', p);
set_db_axis();
ylabel('$|X_{PB}(f)|$ (dB)');
title('Passband spectrum');
tl_labels(tl, '', '$f$ (MHz)');
save_fig(fh, 'q1d');

%% (e) O-QPSK baseband: delay Q by Tsym/2, Q(t) = I(t) for 0 <= t < Tsym/2
d = sps/2;
I = real(xbb);
Q = imag(xbb);
Qo = [I(1:d), Q(1:end-d)];
xoq = I + 1j*Qo;

fh = new_fig('Q1(e)');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_us, real(xoq), 'Color', p.blue);
sym_lines(0:nsym, Tsym_us, p.muted);
set_amp_axis([0 nsym*Tsym_us]);
ylabel('$I(t)$');
nexttile;
plot(t_us, imag(xoq), 'Color', p.orange);
sym_lines(0:nsym, Tsym_us, p.muted);
sym_lines(0.5:nsym, Tsym_us, p.muted, ':');
set_amp_axis([0 nsym*Tsym_us]);
ylabel('$Q(t)$');
tl_labels(tl, 'O-QPSK baseband signal ($Q$ delayed by $T_{sym}/2$)', xl_t);
save_fig(fh, 'q1e');

%% (f) O-QPSK phase jump, 0 <= t < 16 Tsym
% Phase only changes at half-symbol boundaries t = m*Tsym/2
tj_us = (0:2*nsym-1)*Tsym_us/2;
pj_q = phase_jump(xbb, d);
pj_o = phase_jump(xoq, d);

fh = new_fig('Q1(f)', [900 420]);
stem_jump(tj_us, pj_o, p.orange, p.red);
set_phase_axis(nsym*Tsym_us, xl_t);
title('O-QPSK phase jump');
save_fig(fh, 'q1f');

%% (g) QPSK vs O-QPSK phase jump
fh = new_fig('Q1(g)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem_jump(tj_us, pj_q, p.blue, p.red);
set_phase_axis(nsym*Tsym_us, xl_t);
title('QPSK');
nexttile;
stem_jump(tj_us, pj_o, p.orange, p.red);
set_phase_axis(nsym*Tsym_us, xl_t);
title('O-QPSK');
tl_labels(tl, 'Phase jump comparison', '');
save_fig(fh, 'q1g');

fprintf('\nPhase jump count   -pi  -pi/2    0  +pi/2\n');
fprintf('QPSK             %4d  %5d  %3d  %5d\n', jump_count(pj_q));
fprintf('O-QPSK           %4d  %5d  %3d  %5d\n', jump_count(pj_o));

%% Local functions
function sym_lines(k, Tsym_us, color, style)
% Symbol boundaries as vertical lines
if nargin < 4, style = '--'; end
xline(k*Tsym_us, style, 'Color', color, 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function fc_line(f_MHz, label, side, p)
% Carrier marker with its label beside the line
xline(f_MHz, '--', label, 'Color', p.ink, 'LineWidth', 1, ...
    'Interpreter', 'latex', 'FontSize', 15, 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', side);
end

function y = mag_db(X)
y = 20*log10(abs(X)/max(abs(X)) + eps);
end

function pj = phase_jump(x, d)
% Phase jump at every half-symbol boundary, wrapped to [-pi, pi)
% All constellation phases are odd multiples of pi/4, so jumps are multiples of pi/2
idx = 1:d:numel(x);
prev = [x(1), x(idx(2:end)-1)];
k = mod(round((angle(x(idx)) - angle(prev))/(pi/2)) + 2, 4) - 2;
pj = k*pi/2;
end

function c = jump_count(pj)
c = arrayfun(@(v) sum(abs(pj - v) < 1e-9), [-pi -pi/2 0 pi/2]);
end

function stem_jump(t_us, pj, color, hl)
% Stem plot with |jump| = pi highlighted
big = abs(pj) > pi*0.99;
stem(t_us(~big), pj(~big), 'filled', 'Color', color, 'MarkerSize', 5, 'LineWidth', 1.4);
hold on;
if any(big)
    stem(t_us(big), pj(big), 'filled', 'Color', hl, 'MarkerSize', 5, 'LineWidth', 1.4);
    legend('$|\Delta\phi| \leq \pi/2$', '$|\Delta\phi| = \pi$', ...
        'Location', 'northeast', 'Orientation', 'horizontal');
end
end

function set_amp_axis(xl)
xlim(xl); ylim([-1.1 1.1]);
yticks(-1:0.5:1);
end

function set_db_axis()
xlim([-160 160]); ylim([-50 8]);    % 8 dB headroom for labels, range < 60 dB
yticks(-50:10:0);
xticks(-160:40:160);
end

function set_phase_axis(tmax_us, xl)
xlim([0 tmax_us]); ylim([-pi-0.4 pi]);
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'});
xlabel(xl); ylabel('$\Delta\phi$ (rad)');
end
