%% HW3 Q2: random four-path channel and its frequency response
clear; close all; clc;

%% Parameters
ch = hw3_channel();
tau = ch.tau;
Fs = ch.Fs;
ni = ch.ni;
f = linspace(-250e6, 250e6, 10001);    % 50 kHz grid over -250 ~ 250 MHz
f_MHz = f/1e6;

p = hw_palette();
i_idx = 0:3;
xl_f = 'Frequency $f$ (MHz)';

fprintf('m = %d -> rng(%d, ''twister''), gamma = (randn(1,4) + j randn(1,4))/sqrt(2)\n', ch.m, ch.m);

%% (a) Path gains
gm = ch.gamma;
g = ch.g;
h = ch.h;
fprintf('\n(a)  i   %-24s %-24s %-24s\n', 'gamma_i', 'g_i = sqrt(p_i) gamma_i', 'h_i = g_i/sqrt(sum|g|^2)');
for i = 1:4
    fprintf('     %d   %s   %s   %s\n', i_idx(i), cstr(gm(i)), cstr(g(i)), cstr(h(i)));
end
fprintf('\n     i   |gamma_i|   angle(gamma_i)   |g_i|      angle(g_i)   |h_i|      angle(h_i)  (rad)\n');
for i = 1:4
    fprintf('     %d   %.6f   %+.6f        %.6f   %+.6f    %.6f   %+.6f\n', i_idx(i), ...
        abs(gm(i)), angle(gm(i)), abs(g(i)), angle(g(i)), abs(h(i)), angle(h(i)));
end
fprintf('    sum|g_i|^2 = %.8f, sum|h_i|^2 = %.15f\n', sum(abs(g).^2), sum(abs(h).^2));
fprintf('    |h_i|^2 = [%s], p~_i = [%s]\n', sprintf(' %.6f', abs(h).^2), sprintf(' %.6f', ch.pt));
write_q2_table(i_idx, gm, g, h);

%% (b) Continuous frequency response H(f) = sum_i h_i exp(-j2 pi f tau_i)
H = h*exp(-1j*2*pi*tau.'*f);
[Hmax, imax] = max(abs(H));
[Hmin, imin] = min(abs(H));
fprintf('\n(b) |H(f)| over -250 ~ 250 MHz: max = %.4f at %+.3f MHz, min = %.4f (%.1f dB) at %+.3f MHz\n', ...
    Hmax, f_MHz(imax), Hmin, 20*log10(Hmin), f_MHz(imin));
fprintf('    H(0) = %s, |H(0)| = %.6f\n', cstr(H(f == 0)), abs(H(f == 0)));

fh = new_fig('Q2(b)', [900 820]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(f_MHz, abs(H), 'Color', p.blue);
set_f_axis();
ylim([0 2]);
ylabel('$|H(f)|$');
title('Magnitude (linear)');
nexttile;
plot(f_MHz, 20*log10(abs(H)), 'Color', p.blue);
set_f_axis();
ylim([-40 10]);
ylabel('$20\log_{10}|H(f)|$ (dB)');
title('Magnitude (dB)');
nexttile;
plot(f_MHz, angle(H)/pi, 'Color', p.orange, 'LineWidth', 1.2);
set_f_axis();
ylim([-1.1 1.1]);
yticks(-1:0.5:1);
ylabel('$\angle H(f)/\pi$');
title('Phase');
tl_labels(tl, 'Continuous frequency response $H(f) = \sum_i h_i e^{-j2\pi f \tau_i}$', xl_f);
save_fig(fh, 'q2b');

%% (c) Verification and periodicity
% Check 1: tau_i = n_i/Fs, so H(f) equals the DTFT of h[n] at w = 2 pi f/Fs; compare with an FFT
Nfft = 200;    % 1.25 MHz bin spacing
Hfft = fft(ch.hn, Nfft);
mk = -Nfft:Nfft;                          % markers over -250 ~ 250 MHz
f_mk = mk*Fs/Nfft;
H_mk = Hfft(mod(mk, Nfft) + 1);           % FFT bins reused periodically
H_mk_cf = h*exp(-1j*2*pi*tau.'*f_mk);
fprintf('\n(c) check 1: max|FFT_%d{h[n]} - H(f_k)| over %d markers = %.1e\n', Nfft, numel(mk), max(abs(H_mk - H_mk_cf)));

% Check 2: closed-form values
fprintf('    check 2: H(0) - sum(h_i) = %.1e, mean|H(f)|^2 over one period = %.6f (sum|h_i|^2 = 1)\n', ...
    abs(H(f == 0) - sum(h)), mean(abs(H(f >= 0 & f < 250e6)).^2));

% Check 3: periodicity, H(f + P) = H(f) needs P*tau_i integer for every path
tau_ps = round(tau*1e12);
T0 = gcd(gcd(tau_ps(2), tau_ps(3)), tau_ps(4))*1e-12;
fprintf('    check 3: gcd of delays = %g ns -> period 1/gcd = %g MHz, gcd(n_i) = %d\n', ...
    T0*1e9, 1/T0/1e6, gcd(gcd(ni(2), ni(3)), ni(4)));
f0 = linspace(-250e6, 0, 5001);
Hf0 = h*exp(-1j*2*pi*tau.'*f0);
P_try = [250e6 125e6 62.5e6];
dP = zeros(size(P_try));
for k = 1:numel(P_try)
    dP(k) = max(abs(h*exp(-1j*2*pi*tau.'*(f0 + P_try(k))) - Hf0));
    fprintf('             P = %8.4f MHz: max|H(f + P) - H(f)| = %.2e, max||H(f + P)| - |H(f)|| = %.2e\n', ...
        P_try(k)/1e6, dP(k), max(abs(abs(h*exp(-1j*2*pi*tau.'*(f0 + P_try(k)))) - abs(Hf0))));
end
fprintf('    |H(f)| vs |H(-f)|: max difference = %.4f (not symmetric, h_i complex)\n', ...
    max(abs(abs(H) - abs(fliplr(H)))));

fh = new_fig('Q2(c)', [900 820]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(f_MHz, abs(H), 'Color', p.blue); hold on;
plot(f_mk/1e6, abs(H_mk), 'o', 'Color', p.ink, 'MarkerSize', 3.5, 'LineWidth', 0.8);
xline(0, '--', 'Color', p.red, 'LineWidth', 1.1, 'HandleVisibility', 'off');
set_f_axis();
ylim([0 2.4]);
ylabel('$|H(f)|$');
legend('$H(f)$, closed form', sprintf('%d-point FFT of $h[n]$', Nfft), 'Location', 'northeast', 'Orientation', 'horizontal');
title('Magnitude: closed form vs.\ FFT of $h[n]$ (red dashed: $f = 0$, one period on each side)');
nexttile;
plot(f_MHz, angle(H)/pi, 'Color', p.orange, 'LineWidth', 1.2); hold on;
plot(f_mk/1e6, angle(H_mk)/pi, 'o', 'Color', p.ink, 'MarkerSize', 3.5, 'LineWidth', 0.8);
xline(0, '--', 'Color', p.red, 'LineWidth', 1.1);
set_f_axis();
ylim([-1.1 1.1]);
yticks(-1:0.5:1);
ylabel('$\angle H(f)/\pi$');
title('Phase: closed form vs.\ FFT of $h[n]$');
nexttile;
plot(f0/1e6, abs(Hf0), 'Color', p.blue, 'LineWidth', 2.4); hold on;
plot(f0/1e6, abs(h*exp(-1j*2*pi*tau.'*(f0 + 250e6))), '--', 'Color', p.orange, 'LineWidth', 1.4);
set_f_axis();
xlim([-250 0]);
xticks(-250:25:0);
ylim([0 2.4]);
ylabel('$|H|$');
legend('$|H(f)|$', '$|H(f + 250\,\mathrm{MHz})|$', 'Location', 'northeast', 'Orientation', 'horizontal');
title('Left half shifted by $F_s = 250$ MHz lands on the right half');
tl_labels(tl, 'Verification of $H(f)$ and its period $1/\gcd(\tau_i) = 250$ MHz', xl_f);
save_fig(fh, 'q2c');

%% Local functions
function set_f_axis()
xlim([-250 250]);
xticks(-250:50:250);
end

function s = cstr(z)
s = sprintf('%+.6f%+.6fj', real(z), imag(z));
end

function write_q2_table(i_idx, gm, g, h)
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, 'q2a.tex'), 'w');
fprintf(fid, '\\begin{tabular}{c|rr|rr|rr}\n\\toprule\n');
fprintf(fid, '$i$ & \\multicolumn{2}{c|}{$\\gamma_i$} & \\multicolumn{2}{c|}{$g_i = \\sqrt{p_i}\\,\\gamma_i$} & \\multicolumn{2}{c}{$h_i$} \\\\\n');
fprintf(fid, ' & value & $|\\cdot|$, $\\angle\\cdot$ (rad) & value & $|\\cdot|$, $\\angle\\cdot$ (rad) & value & $|\\cdot|$, $\\angle\\cdot$ (rad) \\\\\n\\midrule\n');
for i = 1:4
    fprintf(fid, '%d & %s & %s & %s & %s & %s & %s \\\\\n', i_idx(i), ...
        ctex(gm(i)), pol(gm(i)), ctex(g(i)), pol(g(i)), ctex(h(i)), pol(h(i)));
end
fprintf(fid, '\\midrule\n$\\sum|\\cdot|^2$ & \\multicolumn{2}{c|}{%.6f} & \\multicolumn{2}{c|}{%.6f} & \\multicolumn{2}{c}{%.6f} \\\\\n', ...
    sum(abs(gm).^2), sum(abs(g).^2), sum(abs(h).^2));
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end

function s = ctex(z)
s = sprintf('$%.4f %s %.4fj$', real(z), char(43 + 2*(imag(z) < 0)), abs(imag(z)));
end

function s = pol(z)
s = sprintf('$%.4f$, $%.4f$', abs(z), angle(z));
end
