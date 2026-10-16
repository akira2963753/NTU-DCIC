%% HW2 Q3: m-sequences and Gold code
clear; close all; clc;

%% Parameters
m = 39;    % last two digits of the student ID
L = 63;
r1 = 1 + mod(m, 20);
r2 = 11 + mod(m, 20);

% Taps follow the figure: stage k (left = 1, right = 6) feeds the XOR for each x^k term,
% the XOR result enters stage 1, and stage 6 is the output
taps_top = [1 6];          % p(x)  = x^6 + x + 1                (103_8)
taps_bot = [1 2 5 6];      % p'(x) = x^6 + x^5 + x^2 + x + 1    (147_8)

st1 = dec2bin(r1, 6) - '0';    % MSB -> stage 1
st2 = dec2bin(r2, 6) - '0';
z6 = zeros(1, 6);

p = hw_palette();
l = 0:L-1;

fprintf('m = %d, r1 = %d -> [%s], r2 = %d -> [%s]\n', m, r1, sprintf('%d', st1), r2, sprintf('%d', st2));
fprintf('octal 103 = %s, octal 147 = %s\n', dec2bin(oct2dec(103)), dec2bin(oct2dec(147)));

%% (a) m-sequences: one register loaded, the other all zeros
[s1, top1, bot1] = gold_gen(st1, z6, taps_top, taps_bot, L);
[s2, top2, bot2] = gold_gen(z6, st2, taps_top, taps_bot, L);

fprintf('\n(a) s1: top init = [%s], bottom init = [%s]\n', sprintf('%d', st1), sprintf('%d', z6));
fprintf('    s1[0..62] = %s\n', grp(s1));
fprintf('    s2: top init = [%s], bottom init = [%s]\n', sprintf('%d', z6), sprintf('%d', st2));
fprintf('    s2[0..62] = %s\n', grp(s2));
fprintf('    check: bottom output for s1 all zero = %d, top output for s2 all zero = %d\n', ...
    ~any(bot1), ~any(top2));
inits = {st1, st2};
tps = {taps_top, taps_bot};
seqs = {s1, s2};
for i = 1:2
    [~, states] = lfsr(inits{i}, tps{i}, L + 1);
    nuniq = size(unique(states(1:L, :), 'rows'), 1);
    fprintf('    s%d: ones = %d, zeros = %d, distinct states in 63 clocks = %d, state(63) == state(0): %d\n', ...
        i, sum(seqs{i}), L - sum(seqs{i}), nuniq, isequal(states(L+1, :), inits{i}));
end
write_seq_table(st1, st2, s1, s2);

%% (b) Periodic autocorrelation of d_i = (-1)^s_i
d1 = (-1).^s1;
d2 = (-1).^s2;
R1 = acorr(d1);
R2 = acorr(d2);
fprintf('\n(b) R1[0] = %.4f, R1[l>0] in {%s}\n', R1(1), strjoin(unique(compose('%+.4f', R1(2:end))), ', '));
fprintf('    R2[0] = %.4f, R2[l>0] in {%s}\n', R2(1), strjoin(unique(compose('%+.4f', R2(2:end))), ', '));
fprintf('    -1/63 = %.4f\n', -1/63);

fh = new_fig('Q3(b)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
R = {R1, R2};
for i = 1:2
    nexttile;
    stem(l, R{i}, 'filled', 'Color', p.blue, 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
    yline(-1/L, '--', 'Color', p.red, 'LineWidth', 1);
    set_lag_axis(L); ylim([-0.2 1.15]);
    ylabel(sprintf('$R_%d[l]$', i));
    legend(sprintf('$R_%d[l]$', i), '$-1/63$', 'Location', 'northeast');
    title(sprintf('$d_%d[n] = (-1)^{s_%d[n]}$', i, i));
end
tl_labels(tl, 'Normalized periodic autocorrelation of the m-sequences', 'Lag $l$');
save_fig(fh, 'q3b');

%% (c) Gold code: r1 in the top register and r2 in the bottom register
[s3, top3, bot3] = gold_gen(st1, st2, taps_top, taps_bot, L);
d3 = (-1).^s3;
fprintf('\n(c) s3: top init = [%s], bottom init = [%s]\n', sprintf('%d', st1), sprintf('%d', st2));
fprintf('    s3[0..62] = %s\n', grp(s3));
fprintf('    s3 == s1 xor s2: %d (ones in s3 = %d)\n', isequal(s3, xor(s1, s2)), sum(s3));
write_seq_table3(s3);

C13 = xcorr_p(d1, d3);
C23 = xcorr_p(d2, d3);
C12 = xcorr_p(d1, d2);
tvals = [-17 -1 15];
fprintf('    63*C13 values: %s\n', count_str(C13, L));
fprintf('    63*C23 values: %s\n', count_str(C23, L));
fprintf('    max|C13| = %.4f, max|C23| = %.4f, t(6)/63 = 17/63 = %.4f\n', max(abs(C13)), max(abs(C23)), 17/63);

fh = new_fig('Q3(c)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
C = {C13, C23};
lbl = {'C_{1,3}[l]', 'C_{2,3}[l]'};
for i = 1:2
    nexttile;
    cross_stem(l, C{i}, tvals/L, p);
    set_lag_axis(L); ylim([-0.4 0.4]);
    yticks(tvals/L);
    yticklabels({'$-17/63$', '$-1/63$', '$15/63$'});
    ylabel(['$' lbl{i} '$']);
    title(sprintf('$(u,v) = (%d,3)$', i));
end
tl_labels(tl, 'Periodic cross-correlation with the Gold codeword $d_3[n]$', 'Lag $l$');
save_fig(fh, 'q3c');

%% (d) Check against Gold-code theory
% Preferred pair of degree n = 6 (n mod 4 = 2): cross-correlation takes only
% {-t(n), -1, t(n) - 2}/L with t(n) = 1 + 2^((n+2)/2) = 17
nn = 6;
t = 1 + 2^((nn+2)/2);
fprintf('\n(d) t(6) = %d, allowed values = {%d, -1, %d}/63\n', t, -t, t - 2);
fprintf('    63*C12 values: %s\n', count_str(C12, L));
fprintf('    sum of C over all lags: C12 %.4f, C13 %.4f, C23 %.4f (expected (sum d_u)(sum d_v)/63)\n', ...
    sum(C12), sum(C13), sum(C23));

% Shift-and-add: d1[n] d1[n-l] = d1[n-k] for every l != 0
ok = true;
for ll = 1:L-1
    prod1 = d1.*circshift(d1, ll);
    ok = ok && any(arrayfun(@(k) isequal(prod1, circshift(d1, k)), 0:L-1));
end
fprintf('    shift-and-add holds for d1 at all 62 nonzero lags: %d\n', ok);

% Preferred pair as a decimation: s2 is a cyclic shift of s1[q*n mod 63]
for q = [5 13 21 23]
    dec = s1(mod(q*(0:L-1), L) + 1);
    k = find(arrayfun(@(k) isequal(circshift(dec, k), s2), 0:L-1), 1) - 1;
    if isempty(k), fprintf('    q = %2d: s1[qn] is not a shift of s2\n', q);
    else, fprintf('    q = %2d: s1[qn] = s2 shifted by %d\n', q, k); end
end

% Gold codeword is not an m-sequence: its autocorrelation is three-valued
R3 = acorr(d3);
fprintf('    63*R3 values (l > 0): %s\n', count_str(R3(2:end), L));

% Brute force over every Gold codeword pair
G = zeros(L + 2, L);    % d1, d2 and the 63 codewords d1 * (shifted d2)
G(1, :) = d1;
G(2, :) = d2;
for k = 0:L-1
    G(k+3, :) = (-1).^xor(s1, circshift(s2, -k));
end
vals = [];
for a = 1:size(G, 1)
    for b = a+1:size(G, 1)
        vals = [vals, round(L*xcorr_p(G(a, :), G(b, :)))]; %#ok<AGROW>
    end
end
fprintf('    all %d codeword pairs, all lags: values = {%s}\n', nchoosek(size(G, 1), 2), ...
    sprintf('%d ', unique(vals)));

cnt = [counts(C12, L, tvals); counts(C13, L, tvals); counts(C23, L, tvals)];
e = 2;    % gcd(n, k) for the decimation q = 2^k + 1 = 5
cnt_th = [2^(nn-e-1) - 2^((nn-e)/2-1), 2^nn - 2^(nn-e) - 1, 2^(nn-e-1) + 2^((nn-e)/2-1)];
fprintf('    theoretical counts of {-17, -1, 15} for a preferred pair: %d %d %d\n', cnt_th);

fh = new_fig('Q3(d)', [900 440]);
b = bar(categorical({'-17/63', '-1/63', '15/63'}, {'-17/63', '-1/63', '15/63'}), [cnt; cnt_th].', ...
    'grouped', 'EdgeColor', 'none');
b(1).FaceColor = p.muted; b(2).FaceColor = p.blue; b(3).FaceColor = p.orange; b(4).FaceColor = p.ink;
set(gca, 'TickLabelInterpreter', 'latex');
xticklabels({'$-17/63$', '$-1/63$', '$15/63$'});
for i = 1:4
    xtips = b(i).XEndPoints;
    text(xtips, b(i).YEndPoints + 1.2, string(b(i).YData), 'HorizontalAlignment', 'center', 'FontSize', 11);
end
ylim([0 55]);
ylabel('Number of lags');
legend('$C_{1,2}$', '$C_{1,3}$', '$C_{2,3}$', 'Theory (preferred pair)', 'Location', 'northeast');
title('Cross-correlation values over all 63 lags');
save_fig(fh, 'q3d');

%% Local functions
function [out, states] = lfsr(init, taps, len)
% Fibonacci LFSR: output stage 6, feedback XOR of the tapped stages into stage 1
st = init;
out = zeros(1, len);
states = zeros(len, numel(init));
for n = 1:len
    states(n, :) = st;
    out(n) = st(end);
    fb = mod(sum(st(taps)), 2);
    st = [fb st(1:end-1)];
end
end

function [s, top, bot] = gold_gen(init_top, init_bot, taps_top, taps_bot, len)
top = lfsr(init_top, taps_top, len);
bot = lfsr(init_bot, taps_bot, len);
s = xor(top, bot);
end

function R = acorr(d)
% R[l] = 1/N sum_n d[n] d[(n + l) mod N]
N = numel(d);
R = arrayfun(@(l) sum(d.*circshift(d, -l))/N, 0:N-1);
end

function C = xcorr_p(du, dv)
% C[l] = 1/N sum_n du[n] dv[(n - l) mod N]
N = numel(du);
C = arrayfun(@(l) sum(du.*circshift(dv, l))/N, 0:N-1);
end

function c = counts(C, L, vals)
c = arrayfun(@(v) sum(round(L*C) == v), vals);
end

function s = count_str(C, L)
v = round(L*C);
u = unique(v);
s = strjoin(arrayfun(@(x) sprintf('%d (x%d)', x, sum(v == x)), u, 'UniformOutput', false), ', ');
end

function s = grp(b)
% Bits in groups of 9
s = sprintf('%d', b);
s = strjoin(cellstr(reshape(s, 9, []).'), ' ');
end

function cross_stem(l, C, lv, p)
% Stems colored by their value: -17/63 red, -1/63 gray, 15/63 blue
v = round(63*C);
cols = {p.red, p.muted, p.blue};
for i = 1:3
    idx = v == round(63*lv(i));
    if any(idx)
        stem(l(idx), C(idx), 'filled', 'Color', cols{i}, 'MarkerSize', 4, 'LineWidth', 1.2);
        hold on;
    end
end
yline(lv, ':', 'Color', p.ink2, 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function set_lag_axis(L)
xlim([-1 L]);
xticks(0:6:L-1);
end

function write_seq_table(st1, st2, s1, s2)
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
fid = fopen(fullfile(tabdir, 'q3a.tex'), 'w');
fprintf(fid, '\\begin{tabular}{lcc}\n\\toprule\n');
fprintf(fid, ' & $s_1[n]$ & $s_2[n]$ \\\\\n\\midrule\n');
fprintf(fid, 'Top LFSR initial state & \\texttt{%s} & \\texttt{000000} \\\\\n', sprintf('%d', st1));
fprintf(fid, 'Bottom LFSR initial state & \\texttt{000000} & \\texttt{%s} \\\\\n', sprintf('%d', st2));
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

fid = fopen(fullfile(tabdir, 'q3a_seq.tex'), 'w');
seq_head(fid);
seq_rows(fid, s1, 's_1');
fprintf(fid, '\\midrule\n');
seq_rows(fid, s2, 's_2');
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end

function write_seq_table3(s3)
tabdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'tab');
fid = fopen(fullfile(tabdir, 'q3c_seq.tex'), 'w');
seq_head(fid);
seq_rows(fid, s3, 's_3');
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
end

function seq_head(fid)
fprintf(fid, '\\begin{tabular}{lcc}\n\\toprule\nSequence & $n$ & Bits \\\\\n\\midrule\n');
end

function seq_rows(fid, s, name)
% 63 bits as 7 rows of 9: n = 9r ~ 9r+8
for r = 0:6
    if r == 0, lead = ['$' name '[n]$']; else, lead = ''; end
    fprintf(fid, '%s & %d--%d & \\texttt{%s} \\\\\n', lead, 9*r, 9*r+8, sprintf('%d', s(9*r+(1:9))));
end
end
