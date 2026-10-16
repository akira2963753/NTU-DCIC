%% HW2 Q2: Barker code periodic autocorrelation
clear; close all; clc;

%% Barker codes from Table 1
codes = {
    [-1 -1 -1  1  1 -1  1]
    [-1 -1 -1  1  1  1 -1  1  1 -1  1]
    [-1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1]};
Ns = [7 13];

p = hw_palette();

%% Periodic (and aperiodic, for comparison) autocorrelation
fh = new_fig('Q2', [900 600]);
tl = tiledlayout(numel(Ns), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for N = Ns
    s = codes{cellfun(@numel, codes) == N};
    k = 0:N-1;
    Phi = zeros(1, N);
    for kk = k
        Phi(kk+1) = sum(s.*conj(circshift(s, kk)))/N;    % s*((n - k) mod N)
    end
    A = zeros(1, N);    % aperiodic: 1/N sum_{n=k}^{N-1} s(n) s(n-k)
    for kk = k
        A(kk+1) = sum(s(kk+1:N).*s(1:N-kk))/N;
    end

    fprintf('N = %2d: code = [%s], sum = %+d\n', N, sprintf('%+d ', s), sum(s));
    fprintf('    N*Phi(k), k = 0..%d: %s\n', N-1, sprintf('%+d ', round(N*Phi)));
    fprintf('    off-peak Phi = %s (expected sum^2 - N over N(N-1) = %+.4f)\n', ...
        strjoin(unique(compose('%+.4f', Phi(2:end))), ', '), (sum(s)^2 - N)/(N*(N-1)));
    fprintf('    N*aperiodic(k), k = 0..%d: %s, max off-peak |.| = %d\n', ...
        N-1, sprintf('%+d ', round(N*A)), max(abs(round(N*A(2:end)))));
    write_table(N, round(N*Phi), round(N*A));

    nexttile;
    stem(k, Phi, 'filled', 'Color', p.blue, 'MarkerSize', 5, 'LineWidth', 1.4); hold on;
    plot(k, A, 'o', 'Color', p.orange, 'MarkerSize', 8, 'LineWidth', 1.2);
    yline(Phi(2), '--', 'Color', p.muted, 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim([-0.6 N - 0.4]); xticks(k);
    ylim([-0.4 1.15]);
    yticks(-0.4:0.2:1);
    ylabel(sprintf('$\\Phi_s(k)$, $N = %d$', N));
    legend('Periodic $\Phi_s(k)$', 'Aperiodic (reference)', 'Location', 'northeast');
    title(sprintf('$N = %d$: off-peak $\\Phi_s(k) = %+d/%d$ for all $k \\neq 0$', N, round(N*Phi(2)), N));
end
tl_labels(tl, 'Periodic autocorrelation of the Barker codes', 'Shift $k$');
save_fig(fh, 'q2');

%% Local functions
function write_table(N, P, A)
% N*Phi_s(k) and N*aperiodic(k) as a LaTeX table, needs \usepackage{booktabs}
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, sprintf('q2_n%d.tex', N)), 'w');
fprintf(fid, '\\begin{tabular}{l*{%d}{r}}\n\\toprule\n', N);
fprintf(fid, '$k$%s \\\\\n\\midrule\n', sprintf(' & %d', 0:N-1));
fprintf(fid, '$N\\Phi_s(k)$%s \\\\\n', sprintf(' & $%d$', P));
fprintf(fid, 'Aperiodic $\\times N$%s \\\\\n', sprintf(' & $%d$', A));
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end
