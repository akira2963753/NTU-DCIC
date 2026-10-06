%% HW1 Q3: OTFS baseband
clear; close all; clc;

%% Parameters
seed = 39;
N = 4;    % Doppler / time slots
M = 8;    % delay / subcarriers

p = hw_palette();

%% Delay-Doppler grid: X_DD[l,k] = 1 - 2*b_i, i = k*M + l
b = gen_bits(seed);
Xdd = reshape(1 - 2*b, M, N);    % row l+1, column k+1

fprintf('X_DD[l,k] (rows l = 0~%d, columns k = 0~%d)\n', M-1, N-1);
disp(Xdd);

%% (a) ISFFT: delay-Doppler -> time-frequency
% X_TF[n,m] = 1/sqrt(NM) sum_k sum_l X_DD[l,k] e^{j2pi(nk/N - ml/M)}
n = (0:N-1).';
m = 0:M-1;
l = (0:M-1).';
k = 0:N-1;
Fn = exp(1j*2*pi*n*k/N);        % N x N, [n,k]
Fm = exp(-1j*2*pi*l*m/M);       % M x M, [l,m]
Xtf = Fn*Xdd.'*Fm/sqrt(N*M);    % N x M, [n,m]

% Cross-check with the element-wise double sum
Xtf_sum = zeros(N, M);
for nn = 0:N-1
    for mm = 0:M-1
        e = exp(1j*2*pi*(nn*k/N - mm*l/M));    % M x N, [l,k]
        Xtf_sum(nn+1, mm+1) = sum(Xdd.*e, 'all')/sqrt(N*M);
    end
end
fprintf('(a) max|ISFFT matrix form - double sum| = %.1e\n', max(abs(Xtf - Xtf_sum), [], 'all'));
print_row('X_TF[0,m]', Xtf(1, :));
print_row('X_TF[2,m]', Xtf(3, :));

% Real part on top, imaginary part below; n = 0 left, n = 2 right
parts = {@real, @imag};
names = {'Re', 'Im'};
colors = {p.blue, p.orange};

fh = new_fig('Q3(a)', [900 600]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:2
    ymax = 1.2*max(abs(parts{r}(Xtf([1 3], :))), [], 'all');
    for nn = [0 2]
        nexttile;
        stem_one(m, parts{r}(Xtf(nn+1, :)), colors{r});
        set_idx_axis(M);
        ylim([-ymax ymax]);
        title(['$\mathrm{' names{r} '}\{X_{TF}[' num2str(nn) ',m]\}$']);
    end
end
tl_labels(tl, 'Time-frequency grid after ISFFT', 'Subcarrier index $m$');
save_fig(fh, 'q3a');

%% (b) Unitary IFFT along m -> delay-time, then x[q] with q = n*M + l
% X_DT[n,l] = 1/sqrt(M) sum_m X_TF[n,m] e^{j2pi ml/M}
Xdt = ifft(Xtf, M, 2)*sqrt(M);    % N x M, [n,l]
x = reshape(Xdt.', 1, []);        % q = n*M + l
q = 0:N*M-1;
print_row('x[q]', x);

fh = new_fig('Q3(b)', [900 540]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:2
    nexttile;
    stem_one(q, parts{r}(x), colors{r});
    xline(M*(1:N-1) - 0.5, ':', 'Color', p.muted, 'LineWidth', 1, 'HandleVisibility', 'off');
    set_idx_axis(N*M);
    xticks(0:4:N*M);
    ylim([-2.4 2.4]);
    ylabel(['$\mathrm{' names{r} '}\{x[q]\}$']);
end
tl_labels(tl, 'OTFS time-domain signal (dotted lines separate time slots $n$)', 'Sample index $q = nM + l$');
save_fig(fh, 'q3b');

%% (c) Verification
% 1) Closed form: IFFT along m cancels the FFT along l inside the ISFFT, so
%    X_DT[n,l] = 1/sqrt(N) sum_k X_DD[l,k] e^{j2pi nk/N}  (IDFT along Doppler only)
Xdt_cf = (Fn*Xdd.')/sqrt(N);
x_cf = reshape(Xdt_cf.', 1, []);

% 2) Round trip: receiver FFT along l, then SFFT back to delay-Doppler
Xtf_rx = fft(reshape(x, M, N).', M, 2)/sqrt(M);
Xdd_rx = (Fn'*Xtf_rx*Fm')/sqrt(N*M);
Xdd_rx = Xdd_rx.';
b_rx = (1 - round(real(Xdd_rx(:).')))/2;

fprintf('\n(c) 1) max|x[q] - closed form|       = %.1e\n', max(abs(x - x_cf)));
fprintf('    2) max|X_DD recovered - X_DD|    = %.1e, max|Im| = %.1e, bit errors = %d\n', ...
    max(abs(Xdd_rx - Xdd), [], 'all'), max(abs(imag(Xdd_rx)), [], 'all'), sum(b_rx ~= b));
fprintf('    3) energy: sum|X_DD|^2 = %.4f, sum|X_TF|^2 = %.4f, sum|x|^2 = %.4f\n', ...
    sum(abs(Xdd).^2, 'all'), sum(abs(Xtf).^2, 'all'), sum(abs(x).^2));

fh = new_fig('Q3(c)', [900 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
stem_pair(q, x, p); hold on;
plot(q - 0.12, real(x_cf), 'o', 'Color', p.ink, 'MarkerSize', 8, 'LineWidth', 1);
plot(q + 0.12, imag(x_cf), 'o', 'Color', p.ink, 'MarkerSize', 8, 'LineWidth', 1);
set_idx_axis(N*M);
xticks(0:4:N*M);
ylabel('$x[q]$');
legend('Real (ISFFT + IFFT)', 'Imag (ISFFT + IFFT)', 'Closed form', 'Location', 'eastoutside');
title('Two-step result vs. closed form $\frac{1}{\sqrt{N}}\sum_k X_{DD}[l,k]e^{j2\pi nk/N}$');
nexttile;
i_idx = 0:N*M-1;
stem(i_idx, Xdd(:), 'filled', 'Color', p.blue, 'MarkerSize', 5, 'LineWidth', 1.4); hold on;
plot(i_idx, real(Xdd_rx(:)), 'o', 'Color', p.ink, 'MarkerSize', 8, 'LineWidth', 1);
set_idx_axis(N*M);
xticks(0:4:N*M);
ylim([-1.4 1.4]);
ylabel('$X_{DD}$');
legend('Transmitted', 'Recovered (FFT + SFFT)', 'Location', 'eastoutside');
title('Round trip back to the delay-Doppler grid');
tl_labels(tl, '', 'Index $q = nM + l$ (top), $i = kM + l$ (bottom)');
save_fig(fh, 'q3c');

%% Local functions
function stem_one(idx, v, color)
stem(idx, v, 'filled', 'Color', color, 'MarkerSize', 5, 'LineWidth', 1.4);
end

function stem_pair(idx, z, p)
% Real and imaginary parts as side-by-side stems
stem(idx - 0.12, real(z), 'filled', 'Color', p.blue, 'MarkerSize', 5, 'LineWidth', 1.4);
hold on;
stem(idx + 0.12, imag(z), 'filled', 'Color', p.orange, 'MarkerSize', 5, 'LineWidth', 1.4);
legend('Real', 'Imag', 'Location', 'eastoutside');
end

function set_idx_axis(len)
xlim([-0.6 len - 0.4]);
xticks(0:len-1);
ylim('padded');
end

function print_row(name, z)
fprintf('%s =\n', name);
fprintf('    %+.4f%+.4fj\n', [real(z); imag(z)]);
end
