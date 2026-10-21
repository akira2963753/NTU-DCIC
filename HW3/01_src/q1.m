%% HW3 Q1: power delay profile statistics
clear; close all; clc;

%% Parameters
ch = hw3_channel();
P_dB = ch.P_dB;
tau = ch.tau;
tau_ns = tau*1e9;
Fs = ch.Fs;
Bs = [2.5e6 250e6];    % signal bandwidths of Q1(c)

p = hw_palette();
i_idx = 0:3;

%% (a) Linear and normalized powers
pl = ch.p;
pt = ch.pt;
fprintf('(a) p   = [%s]\n', sprintf(' %.8f', pl));
fprintf('    p~  = [%s]\n', sprintf(' %.8f', pt));
fprintf('    sum(p) = %.8g, sum(p~) = %.8g\n', sum(pl), sum(pt));

%% (b) Delay moments
m1_i = pt.*tau_ns;              % contribution to the first moment (ns)
m2_i = pt.*tau_ns.^2;           % contribution to the second raw moment (ns^2)
tau_bar = sum(m1_i);
c2_i = pt.*(tau_ns - tau_bar).^2;    % contribution to the second central moment (ns^2)
tau_rms = sqrt(sum(c2_i));
tau_max = max(tau_ns) - min(tau_ns);

fprintf('\n(b)  i   tau(ns)   p~_i       p~_i*tau_i(ns)  p~_i*tau_i^2(ns^2)  p~_i*(tau_i-mean)^2(ns^2)\n');
for i = 1:4
    fprintf('     %d   %5.0f   %.6f   %12.6f   %16.6f   %20.6f\n', i_idx(i), tau_ns(i), pt(i), m1_i(i), m2_i(i), c2_i(i));
end
fprintf('    sum                    %12.6f   %16.6f   %20.6f\n', sum(m1_i), sum(m2_i), sum(c2_i));
fprintf('    mean excess delay = %.6f ns\n', tau_bar);
fprintf('    RMS delay spread  = %.6f ns (check: sqrt(E[tau^2] - mean^2) = %.6f ns)\n', ...
    tau_rms, sqrt(sum(m2_i) - tau_bar^2));
fprintf('    max excess delay  = %.0f ns\n', tau_max);
write_q1_table(i_idx, P_dB, pl, pt, tau_ns, m1_i, m2_i, c2_i);

fh = new_fig('Q1(b)', [900 720]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1 2]);
xregion(tau_bar - tau_rms, tau_bar + tau_rms, 'FaceColor', p.aqua, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
hold on;
stem(tau_ns, pt, 'filled', 'Color', p.blue, 'MarkerSize', 6, 'LineWidth', 2);
xline(tau_bar, '--', 'Color', p.red, 'LineWidth', 1.3);
text(tau_ns, pt + 0.03, compose('$%.4f$', pt), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(tau_bar + 3, 0.55, sprintf('$\\bar\\tau = %.2f$ ns', tau_bar), 'Color', p.red, 'FontSize', 12);
text(tau_bar + tau_rms + 3, 0.43, sprintf('$\\bar\\tau \\pm \\tau_{\\mathrm{RMS}}$, $\\tau_{\\mathrm{RMS}} = %.2f$ ns', tau_rms), ...
    'Color', p.aqua, 'FontSize', 12);
xlim([-30 230]); ylim([0 0.6]);
xticks(0:20:220);
xlabel('Excess delay $\tau$ (ns)');
ylabel('$\tilde p_i$');
title('Normalized power delay profile');
nexttile;
bar(i_idx, m1_i, 0.55, 'FaceColor', p.blue, 'EdgeColor', 'none');
text(i_idx, m1_i, compose('$%.2f$', m1_i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
ylim([0 16]);
xlabel('Path $i$');
ylabel('$\tilde p_i \tau_i$ (ns)');
title(sprintf('First moment ($\\bar\\tau = %.2f$ ns)', tau_bar));
nexttile;
bar(i_idx, c2_i, 0.55, 'FaceColor', p.orange, 'EdgeColor', 'none');
text(i_idx, c2_i, compose('$%.0f$', c2_i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
ylim([0 2000]);
xlabel('Path $i$');
ylabel('$\tilde p_i (\tau_i - \bar\tau)^2$ (ns$^2$)');
title(sprintf('Second central moment ($\\tau_{\\mathrm{RMS}}^2 = %.0f$ ns$^2$)', tau_rms^2));
tl_labels(tl, 'Four-path channel: power delay profile and delay moments', '');
save_fig(fh, 'q1b');

%% (c) Coherence bandwidth
Bc = 1/(5*tau_rms*1e-9);
fprintf('\n(c) Bc = 1/(5 tau_rms) = %.6f MHz\n', Bc/1e6);
for B = Bs
    if B < Bc, cls = 'flat fading'; else, cls = 'frequency selective'; end
    fprintf('    Bs = %6.1f MHz: Bs/Bc = %8.4f -> %s\n', B/1e6, B/Bc, cls);
end
fprintf('    (K = 50 in Bc = 1/(K tau_rms) would give Bc = %.4f MHz)\n', 1/(50*tau_rms*1e-9)/1e6);

%% (d) Integer sample delays at Fs = 250 MHz
ni_raw = tau*Fs;
ni = ch.ni;
L = ch.L;
fprintf('\n(d) Ts = 1/Fs = %g ns\n', 1e9/Fs);
fprintf('    tau_i*Fs = [%s] -> n_i = [%s], max rounding error = %.1e\n', ...
    sprintf(' %.4f', ni_raw), sprintf(' %d', ni), max(abs(ni_raw - ni)));
fprintf('    discrete channel length L = max(n_i) + 1 = %d taps (n = 0 ~ %d)\n', L, L - 1);
fprintf('    minimum CP length = L - 1 = %d samples = %g ns\n', L - 1, (L - 1)/Fs*1e9);

fh = new_fig('Q1(d)', [900 400]);
n = 0:L-1;
pn = zeros(1, L);
pn(ni + 1) = pt;
stem(n, pn, 'filled', 'Color', p.blue, 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
stem(ni, pt, 'filled', 'Color', p.red, 'MarkerSize', 6, 'LineWidth', 2);
for i = 1:4
    text(ni(i), pt(i) + 0.04, sprintf('$n_%d = %d$', i_idx(i), ni(i)), 'Color', p.red, 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end
xlim([-2 L + 1]); ylim([0 0.6]);
xticks(0:5:L-1);
xlabel('Sample index $n$ ($T_s = 4$ ns)');
ylabel('$\tilde p$ at tap $n$');
title(sprintf('Sampled power delay profile: $L = %d$ taps, minimum CP length $L - 1 = %d$ samples', L, L - 1));
save_fig(fh, 'q1d');

%% Local functions
function write_q1_table(i_idx, P_dB, pl, pt, tau_ns, m1_i, m2_i, c2_i)
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, 'q1.tex'), 'w');
fprintf(fid, '\\begin{tabular}{crcccrrr}\n\\toprule\n');
fprintf(fid, ['$i$ & $P_{i,\\mathrm{dB}}$ & $p_i$ & $\\tilde p_i$ & $\\tau_i$ (ns) & ' ...
    '$\\tilde p_i \\tau_i$ (ns) & $\\tilde p_i \\tau_i^2$ (ns$^2$) & $\\tilde p_i (\\tau_i - \\bar\\tau)^2$ (ns$^2$) \\\\\n\\midrule\n']);
for i = 1:4
    fprintf(fid, '%d & $%d$ dB & %s & %s & %d & %.4f & %.2f & %.2f \\\\\n', i_idx(i), P_dB(i), ...
        sig6(pl(i)), sig6(pt(i)), round(tau_ns(i)), m1_i(i), m2_i(i), c2_i(i));
end
fprintf(fid, '\\midrule\nSum & & %s & %s & & %.4f & %.2f & %.2f \\\\\n', ...
    sig6(sum(pl)), sig6(sum(pt)), sum(m1_i), sum(m2_i), sum(c2_i));
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end

function s = sig6(x)
% Six significant digits, keeping trailing zeros
s = sprintf('%.*f', max(0, 5 - floor(log10(abs(x)))), x);
end
