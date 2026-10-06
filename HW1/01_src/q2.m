%% HW1 Q2: OFDM baseband
clear; close all; clc;

%% Parameters
seed = 39;
N = 16;
Ts = 25e-9;
os = 16;    % t = u*Ts/os, u = 0 ~ os*N

p = hw_palette();
xl_t = '$t$ (ns)';

%% (a) Useful symbol duration and subcarrier spacing
Tsym = N*Ts;
df = 1/Tsym;
fprintf('(a) Tsym = N*Ts = %g ns, df = 1/Tsym = %g MHz\n', Tsym*1e9, df/1e6);

%% (b) QPSK mapping (white points in Fig. 1) and subcarrier table
b = gen_bits(seed);
b0 = b(1:2:end);
b1 = b(2:2:end);

% (0,0)->270, (0,1)->0, (1,0)->180, (1,1)->90 deg
lut = [-1j, 1, -1, 1j];
X = lut(2*b0 + b1 + 1);

k = 0:N-1;
fk = k*df;
fk(k >= N/2) = (k(k >= N/2) - N)*df;

fprintf('\n(b)  k  (b2k,b2k+1)  X_k   f_k (MHz)\n');
for i = 1:N
    fprintf('    %2d     (%d,%d)     %-3s   %+6.1f\n', k(i), b0(i), b1(i), cplx_str(X(i)), fk(i)/1e6);
end
write_table(k, b0, b1, X, fk);

%% (c) Individual subcarriers k = 1, 14
u = 0:os*N;
t = u*Ts/os;
t_ns = t*1e9;
ks = [1 14];

fh = new_fig('Q2(c)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
fprintf('\n(c) Cycle check over 0 <= t <= Tsym\n');
for kk = ks
    yk = X(kk+1)*exp(1j*2*pi*fk(kk+1)*t);

    % Expected cycles = f_k*Tsym; measured from the total unwrapped phase rotation
    cyc_exp = fk(kk+1)*Tsym;
    ph = unwrap(angle(yk));
    cyc_meas = (ph(end) - ph(1))/(2*pi);
    r = real(yk);
    r = r(abs(r) > 1e-9);
    nzc = sum(diff(sign(r)) ~= 0);
    fprintf('    k = %2d: f_k = %+g MHz, expected %+g, measured %+.4f cycles, ', ...
        kk, fk(kk+1)/1e6, cyc_exp, cyc_meas);
    fprintf('Re zero crossings = %d (expect %d), |y(Tsym) - y(0)| = %.1e\n', ...
        nzc, 2*abs(cyc_exp), abs(yk(end) - yk(1)));

    nexttile;
    plot(t_ns, real(yk), 'Color', p.blue); hold on;
    plot(t_ns, imag(yk), '--', 'Color', p.orange);
    xline((0:abs(cyc_exp))*Tsym/abs(cyc_exp)*1e9, ':', 'Color', p.muted, ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    set_amp_axis([0 Tsym*1e9]);
    ylabel(sprintf('$y_{%d}(t)$', kk));
    legend('Real', 'Imag', 'Location', 'eastoutside');
    title(sprintf('$k = %d$: $f_k = %+g$ MHz, $X_k = %s$, %s in $T_{sym}$', ...
        kk, fk(kk+1)/1e6, cplx_str(X(kk+1)), plural(abs(cyc_exp), 'cycle')));
end
tl_labels(tl, '', xl_t);
save_fig(fh, 'q2c');

%% (d) Sum of all 16 subcarriers
y = X*exp(1j*2*pi*fk.'*t);

fh = new_fig('Q2(d)');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_ns, real(y), 'Color', p.blue);
set_sum_axis(Tsym);
ylabel('$\mathrm{Re}\{y(t)\}$');
nexttile;
plot(t_ns, imag(y), 'Color', p.orange);
set_sum_axis(Tsym);
ylabel('$\mathrm{Im}\{y(t)\}$');
tl_labels(tl, 'OFDM symbol $y(t) = \sum_{k=0}^{N-1} X_k e^{j2\pi f_k t}$', xl_t);
save_fig(fh, 'q2d');

%% (e) 16-point IFFT
x = ifft(X);    % ifft already includes the 1/N factor

% Cross-check against the definition and against y(nTs)
n = (0:N-1).';
x_def = (exp(1j*2*pi*n*k/N)*X.').'/N;
y_n = y(1:os:os*N);
fprintf('\n(e) max|ifft - definition| = %.1e, max|N*x[n] - y(nTs)| = %.1e\n', ...
    max(abs(x - x_def)), max(abs(N*x - y_n)));

tn_ns = (0:N-1)*Ts*1e9;
fh = new_fig('Q2(e)', [900 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_ns, real(y), 'Color', p.blue); hold on;
stem(tn_ns, real(N*x), 'filled', 'Color', p.ink, 'MarkerSize', 5, 'LineWidth', 1.2);
set_sum_axis(Tsym);
ylabel('Real part');
legend('$\mathrm{Re}\{y(t)\}$', '$\mathrm{Re}\{Nx[n]\}$', 'Location', 'eastoutside');
nexttile;
plot(t_ns, imag(y), 'Color', p.orange); hold on;
stem(tn_ns, imag(N*x), 'filled', 'Color', p.ink, 'MarkerSize', 5, 'LineWidth', 1.2);
set_sum_axis(Tsym);
ylabel('Imaginary part');
legend('$\mathrm{Im}\{y(t)\}$', '$\mathrm{Im}\{Nx[n]\}$', 'Location', 'eastoutside');
tl_labels(tl, '16-point IFFT samples $Nx[n]$ at $t = nT_s$ vs. $y(t)$', xl_t);
save_fig(fh, 'q2e');

%% (f) Signed f_k = (k-N)df vs. unsigned f_k = k*df
% At t = nTs: e^{j2pi(k-N)df nTs} = e^{j2pi kn/N} e^{-j2pi n} = e^{j2pi kn/N}
y_uns = X*exp(1j*2*pi*(k*df).'*t);
fprintf('(f) max|y_signed(nTs) - y_unsigned(nTs)| = %.1e, ', max(abs(y_n - y_uns(1:os:os*N))));
fprintf('max|y_signed(t) - y_unsigned(t)| over all t = %.2f\n', max(abs(y - y_uns)));

fh = new_fig('Q2(f)', [900 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_ns, real(y), 'Color', p.blue); hold on;
plot(t_ns, real(y_uns), '--', 'Color', p.muted, 'LineWidth', 1.2);
plot(tn_ns, real(N*x), 'o', 'Color', p.ink, 'MarkerFaceColor', p.ink, 'MarkerSize', 5);
set_sum_axis(Tsym);
ylabel('Real part');
legend('$f_k = (k-N)\Delta f$', '$f_k = k\Delta f$', '$Nx[n]$', 'Location', 'eastoutside');
nexttile;
plot(t_ns, imag(y), 'Color', p.orange); hold on;
plot(t_ns, imag(y_uns), '--', 'Color', p.muted, 'LineWidth', 1.2);
plot(tn_ns, imag(N*x), 'o', 'Color', p.ink, 'MarkerFaceColor', p.ink, 'MarkerSize', 5);
set_sum_axis(Tsym);
ylabel('Imaginary part');
legend('$f_k = (k-N)\Delta f$', '$f_k = k\Delta f$', '$Nx[n]$', 'Location', 'eastoutside');
tl_labels(tl, 'Signed vs. unsigned subcarrier frequencies agree only at $t = nT_s$', xl_t);
save_fig(fh, 'q2f');

%% Local functions
function s = cplx_str(z)
% X_k in {1, j, -1, -j} as text
names = {'1', 'j', '-1', '-j'};
s = names{mod(round(angle(z)/(pi/2)), 4) + 1};
end

function write_table(k, b0, b1, X, fk)
% LaTeX table for Q2(b) as two side-by-side halves, needs \usepackage{booktabs}
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, 'q2b.tex'), 'w');
half = numel(k)/2;
for h = 0:1
    if h == 1, fprintf(fid, '\\hspace{2em}\n'); end
    fprintf(fid, '\\begin{tabular}{cccr}\n\\toprule\n');
    fprintf(fid, '$k$ & $(b_{2k}, b_{2k+1})$ & $X_k$ & $f_k$ (MHz) \\\\\n\\midrule\n');
    for i = h*half + (1:half)
        fprintf(fid, '%d & (%d, %d) & $%s$ & $%s$ \\\\\n', k(i), b0(i), b1(i), cplx_str(X(i)), ...
            strrep(sprintf('%+.1f', fk(i)/1e6), '+0.0', '0.0'));
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
end
fclose(fid);
end

function set_amp_axis(xl)
xlim(xl); ylim([-1.2 1.2]);
yticks(-1:0.5:1);
end

function s = plural(n, word)
s = sprintf('%d %s', n, word);
if n ~= 1, s = [s 's']; end
end

function set_sum_axis(Tsym)
xlim([0 Tsym*1e9]);
xticks(0:50:Tsym*1e9);
ylim('padded');
end
