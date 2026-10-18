%% HW2 Q1: Zadoff-Chu sequence
clear; close all; clc;

%% Parameters
m = 39;    % last two digits of the student ID
N = 53;
u = 1 + mod(m, 52);
v = 1 + mod(u + 16, 52);
u0 = 3 + mod(7*m + 2, 48);
tau0 = mod(5*m + 3, 48);
alpha = exp(1j*pi/5);
Omega = [u0 - 2, u0, u0 + 2];

p = hw_palette();
n = 0:N-1;
l = 0:N-1;
xl_n = 'Sample index $n$';
xl_l = 'Lag $l$';

fprintf('m = %d, N = %d, u = %d, v = %d, u0 = %d, tau0 = %d, Omega = {%d, %d, %d}\n', ...
    m, N, u, v, u0, tau0, Omega);
fprintf('gcd(u,N) = %d, gcd(v,N) = %d, gcd(u-v,N) = %d\n', gcd(u, N), gcd(v, N), gcd(abs(u - v), N));

%% (a) Root-u Zadoff-Chu sequence
xu = zc(u, N);
ph = unwrap(angle(xu));
ph_true = -pi*u*n.*(n+1)/N;    % analytic (continuous) phase

fprintf('\n(a) max||x_u| - 1| = %.1e\n', max(abs(abs(xu) - 1)));
fprintf('    max|wrap(unwrap phase - analytic phase)| = %.1e\n', max(abs(wrap_pi(ph - ph_true))));
fprintf('    unwrap phase range = [%.2f, %.2f] rad, analytic phase range = [%.2f, %.2f] rad\n', ...
    min(ph), max(ph), min(ph_true), max(ph_true));
fprintf('    first samples x_u[0..4]:\n');
fprintf('    %+.4f%+.4fj\n', [real(xu(1:5)); imag(xu(1:5))]);

fh = new_fig('Q1(a)', [900 900]);
tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem_one(n, real(xu), p.blue);
set_idx_axis(N, 8); ylim([-1.2 1.2]);
title('$\mathrm{Re}\{x_u[n]\}$');
nexttile;
stem_one(n, imag(xu), p.orange);
set_idx_axis(N, 8); ylim([-1.2 1.2]);
title('$\mathrm{Im}\{x_u[n]\}$');
nexttile;
stem_one(n, abs(xu), p.aqua);
set_idx_axis(N, 8); ylim([0 1.2]);
title('$|x_u[n]|$');
nexttile;
plot(n, ph, '.-', 'Color', p.ink2, 'MarkerSize', 12, 'LineWidth', 1.1);
set_idx_axis(N, 8);
ylim('padded');
title('Unwrapped phase $\angle x_u[n]$ (rad)');
nexttile([1 2]);
% Shift each unwrap sample by its own 2*pi*k: lands exactly on the analytic parabola
ph_fit = ph + 2*pi*round((ph_true - ph)/(2*pi));
plot(n, ph_true, '-', 'Color', p.ink2, 'LineWidth', 1.2); hold on;
plot(n, ph_fit, 'o', 'Color', p.orange, 'MarkerSize', 6, 'LineWidth', 1.1);
set_idx_axis(N, 4);
ylim('padded');
ylabel('rad');
legend('Analytic $-\pi u n(n+1)/N$', 'Unwrapped phase $+\,2\pi k_n$', 'Location', 'southwest');
title('Analytic (continuous) phase');
tl_labels(tl, sprintf('Root-%d Zadoff-Chu sequence ($N = %d$)', u, N), xl_n);
save_fig(fh, 'q1a');

%% (b) Periodic vs. linear autocorrelation
Ru = pcorr(xu, xu);
Rlin = lcorr(xu);

% Closed form of the linear version: the summand is c_l e^{-j2pi u l n/N},
% so the partial sum over n = l ~ N-1 is a truncated geometric series
Rlin_cf = abs(sin(pi*u*l.*(N-l)/N)./sin(pi*u*l/N))/N;
Rlin_cf(1) = 1;

fprintf('\n(b) R_u[0] = %.4f%+.4fj, max_{l>0}|R_u[l]| = %.1e\n', real(Ru(1)), imag(Ru(1)), max(abs(Ru(2:end))));
[pk, ipk] = max(abs(Rlin(2:end)));
fprintf('    linear: max_{l>0}|R_lin[l]| = %.4f at l = %d, mean_{l>0}|R_lin| = %.4f\n', ...
    pk, ipk, mean(abs(Rlin(2:end))));
fprintf('    max|(|R_lin|) - closed form| = %.1e\n', max(abs(abs(Rlin) - Rlin_cf)));
fprintf('    |R_lin[l]| for l = 1..6: %s\n', sprintf('%.4f ', abs(Rlin(2:7))));

fh = new_fig('Q1(b)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem_one(l, abs(Ru), p.blue);
set_idx_axis(N); ylim([0 1.15]);
ylabel('$|R_u[l]|$');
title('Periodic autocorrelation');
nexttile;
stem_one(l, abs(Rlin), p.orange); hold on;
plot(l, Rlin_cf, 'o', 'Color', p.ink, 'MarkerSize', 7, 'LineWidth', 0.9);
set_idx_axis(N); ylim([0 1.15]);
ylabel('$|R_u^{\mathrm{lin}}[l]|$');
legend('Simulated', 'Closed form', 'Location', 'northeast');
title('Linear (non-periodic) autocorrelation');
tl_labels(tl, sprintf('Normalized autocorrelation of $x_{%d}[n]$', u), xl_l);
save_fig(fh, 'q1b');

%% (c) Periodic cross-correlation between roots u and v
xv = zc(v, N);
Cuv = pcorr(xu, xv);
fprintf('\n(c) |C_uv[l]|: min = %.6f, max = %.6f, 1/sqrt(N) = %.6f\n', ...
    min(abs(Cuv)), max(abs(Cuv)), 1/sqrt(N));

fh = new_fig('Q1(c)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem_one(l, abs(Cuv), p.blue);
yline(1/sqrt(N), '--', 'Color', p.red, 'LineWidth', 1.1);
legend(sprintf('$|C_{%d,%d}[l]|$', u, v), '$1/\sqrt{N} = 0.1374$', 'Location', 'northeast', 'Orientation', 'horizontal');
set_idx_axis(N); ylim([0 0.2]);
ylabel(sprintf('$|C_{%d,%d}[l]|$', u, v));
title('Magnitude');
nexttile;
stem_one(l, angle(Cuv)/pi, p.orange);
set_idx_axis(N); ylim([-1.1 1.1]);
yticks(-1:0.5:1);
ylabel(sprintf('$\\angle C_{%d,%d}[l]/\\pi$', u, v));
title('Phase');
tl_labels(tl, sprintf('Periodic cross-correlation between roots $u = %d$ and $v = %d$', u, v), xl_l);
save_fig(fh, 'q1c');

%% (d) Root and delay detection
r = alpha*circshift(zc(u0, N), tau0);    % r[n] = alpha x_u0[(n - tau0) mod N]
tau = 0:N-1;
D = zeros(numel(Omega), N);
Cpk = zeros(1, numel(Omega));
for i = 1:numel(Omega)
    c = pcorr(r, zc(Omega(i), N));
    D(i, :) = abs(c);
    [~, k] = max(D(i, :));
    Cpk(i) = c(k);
end
[Dmax, idx] = max(D(:));
[iu, it] = ind2sub(size(D), idx);
u_hat = Omega(iu);
tau_hat = tau(it);

fprintf('\n(d) r[n] = exp(j pi/5) x_%d[(n - %d) mod %d]\n', u0, tau0, N);
for i = 1:numel(Omega)
    [pk, k] = max(D(i, :));
    fprintf('    u = %2d: max D = %.4f at tau = %2d, min D = %.4f, ', Omega(i), pk, tau(k), min(D(i, :)));
    fprintf('mean D = %.4f\n', mean(D(i, :)));
end
fprintf('    detected (u, tau) = (%d, %d), true (u0, tau0) = (%d, %d) -> match = %d\n', ...
    u_hat, tau_hat, u0, tau0, u_hat == u0 && tau_hat == tau0);
fprintf('    peak D = %.4f, phase of the peak correlation = %.4f*pi (alpha = exp(j pi/5) -> 0.2*pi)\n', ...
    Dmax, angle(Cpk(iu))/pi);
D2 = sort(D(:), 'descend');
fprintf('    second largest D over all (u, tau) = %.4f, peak-to-second ratio = %.2f (%.1f dB)\n', ...
    D2(2), D2(1)/D2(2), 20*log10(D2(1)/D2(2)));

fh = new_fig('Q1(d)', [900 760]);
tl = tiledlayout(numel(Omega), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(Omega)
    nexttile;
    stem_one(tau, D(i, :), p.blue); hold on;
    if i == iu
        stem(tau_hat, D(i, it), 'filled', 'Color', p.red, 'MarkerSize', 6, 'LineWidth', 1.8);
        text(tau_hat + 1, 0.9, sprintf('peak $D = %.3f$ at $\\tau = %d$', D(i, it), tau_hat), ...
            'Color', p.red, 'FontSize', 13);
    end
    yline(1/sqrt(N), '--', 'Color', p.muted, 'LineWidth', 1, 'HandleVisibility', 'off');
    set_idx_axis(N); ylim([0 1.15]);
    ylabel(sprintf('$D(%d,\\tau)$', Omega(i)));
    if Omega(i) == u0
        title(sprintf('$u = %d = u_0$', Omega(i)));
    else
        title(sprintf('$u = %d$', Omega(i)));
    end
end
tl_labels(tl, sprintf('Detection metric $D(u,\\tau)$ (dashed: $1/\\sqrt{N}$)'), 'Cyclic delay $\tau$');
save_fig(fh, 'q1d');

%% Local functions
function x = zc(u, N)
% Root-u Zadoff-Chu sequence of odd length N
n = 0:N-1;
x = exp(-1j*pi*u*n.*(n+1)/N);
end

function R = pcorr(a, b)
% R[l] = 1/N sum_n a[n] b*[(n - l) mod N], l = 0 ~ N-1
N = numel(a);
R = zeros(1, N);
for l = 0:N-1
    R(l+1) = sum(a.*conj(circshift(b, l)))/N;
end
end

function R = lcorr(a)
% R_lin[l] = 1/N sum_{n=l}^{N-1} a[n] a*[n - l], l = 0 ~ N-1
N = numel(a);
R = zeros(1, N);
for l = 0:N-1
    n = l:N-1;
    R(l+1) = sum(a(n+1).*conj(a(n-l+1)))/N;
end
end

function y = wrap_pi(x)
y = mod(x + pi, 2*pi) - pi;
end

function stem_one(idx, v, color)
stem(idx, v, 'filled', 'Color', color, 'MarkerSize', 4, 'LineWidth', 1.2);
end

function set_idx_axis(len, step)
if nargin < 2, step = 4; end
xlim([-1 len]);
xticks(0:step:len-1);
end
