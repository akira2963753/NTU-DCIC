%% HW4 Q3: hardware symbol boundary detector (threshold, word lengths, test vectors, comparison)
% Phi_NM(m) = |sum z(m-r) z*(m-r-L)|^2 / (sum |z(m-r)|^2)^2 > eta  (corrected Q3 metric)
clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'plot_util'));

%% Parameters
d = hw4_data();
P = sbd_param();
R = P.R; L = P.L; Ns = d.Ns;
p = hw_palette();

z = d.z;                                     % preamble discarded, as in Q2(c)
zn = d.zn;                                   % the same noisy record as Q2(d)
M = numel(z);
nsym = M/Ns;
bnd = (1:nsym)*Ns;                           % symbol ends 160, 320, 480
far = true(M, 1);                            % samples outside every CP correlation lobe
for b = bnd, far(abs((1:M) - b) < R) = false; end
[~, nm_c] = td_metric(z, R, L);
[~, nm_n] = td_metric(zn, R, L);
cases = {'clean', z, nm_c, 'noiseless'; 'noisy', zn, nm_n, sprintf('SNR = %d dB', d.snr_db)};

%% (a) Threshold eta
ntrial = 2000;
rng(2);
pk_mc = zeros(ntrial, nsym);                 % largest Phi_NM inside each lobe
far_mc = zeros(ntrial, 1);                   % largest Phi_NM outside the lobes
ndet_mc = zeros(ntrial, 1);
for t = 1:ntrial
    [~, nm_t] = td_metric(z + d.noise(), R, L);
    pk_mc(t, :) = lobe_peaks(nm_t, bnd, R, M);
    far_mc(t) = max(nm_t(far));
    ndet_mc(t) = numel(sbd_peak_search(nm_t, P.ETA, P.W_WIN));
end
fprintf('(a) R = %d, L = %d, lobe = |m - symbol end| < R\n', R, L);
fprintf('    noiseless : lobe peaks = %s, max outside lobes = %.4f (at m = %d)\n', ...
    sprintf(' %.4f', lobe_peaks(nm_c, bnd, R, M)), max(nm_c(far)), find(nm_c == max(nm_c(far)) & far, 1));
fprintf('    SNR 12 dB : lobe peaks = %s, max outside lobes = %.4f\n', ...
    sprintf(' %.4f', lobe_peaks(nm_n, bnd, R, M)), max(nm_n(far)));
fprintf('    Monte Carlo (%d): min lobe peak = %s, 0.1%% quantile = %s\n', ntrial, ...
    sprintf(' %.4f', min(pk_mc)), sprintf(' %.4f', quantile(pk_mc, 0.001)));
fprintf('                      max outside lobes = %.4f, 99.9%% quantile = %.4f\n', max(far_mc), quantile(far_mc, 0.999));
fprintf('    -> eta = %.4f (U2.14 code %d = 0x%04X), margin to lobe peak %.4f, to outside max %.4f\n', ...
    P.ETA, P.ETA_Q, P.ETA_Q, min(pk_mc(:)) - P.ETA, P.ETA - max(far_mc));
fprintf('    detections per record with eta = %.2f, W = %d (float): %s\n', P.ETA, P.W_WIN, ...
    sprintf(' %d:%d', [unique(ndet_mc) histc(ndet_mc, unique(ndet_mc))].'));

fh = new_fig('Q3(a)', [900 760]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
for b = bnd, xregion(b - R + 1, min(M, b + R - 1), 'FaceColor', p.aqua, 'FaceAlpha', 0.10, 'EdgeColor', 'none'); hold on; end
h1 = plot(1:M, nm_c, ':', 'Color', p.muted, 'LineWidth', 1.2);
h2 = plot(1:M, nm_n, 'Color', p.blue);
h3 = yline(P.ETA, '--', 'Color', p.red, 'LineWidth', 1.3, 'Layer', 'bottom');
xlim([1 M]); ylim([0 1.25]); xticks([1 40:40:M]);
legend([h1 h2 h3], {'noiseless', sprintf('SNR = %d dB', d.snr_db), sprintf('$\\eta = %.2f$', P.ETA)}, ...
    'Location', 'north', 'Orientation', 'horizontal');
xlabel('$m$'); ylabel('$\Phi_{NM}(m)$');
title('$\Phi_{NM}(m)$ and the threshold (shaded: CP correlation lobes, $|m - 160k| < R$)');
nexttile;
edges = 0:0.02:1.3;
histogram(far_mc, edges, 'Normalization', 'probability', 'FaceColor', p.orange, 'FaceAlpha', 0.6, 'EdgeColor', 'none'); hold on;
histogram(pk_mc(:), edges, 'Normalization', 'probability', 'FaceColor', p.blue, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xline(P.ETA, '--', 'Color', p.red, 'LineWidth', 1.3);
xlim([0 1.3]);
legend({'max outside the lobes', 'peak inside each lobe', sprintf('$\\eta = %.2f$', P.ETA)}, 'Location', 'northeast');
xlabel('$\Phi_{NM}$'); ylabel('Probability');
title(sprintf('%d noise realizations at SNR = %d dB', ntrial, d.snr_db));
save_fig(fh, 'q3a');

%% (c) Input quantization and word lengths
amax = max(abs([real(z); imag(z); real(zn); imag(zn)]));
fprintf('\n(c) max |Re|,|Im| of the input (noiseless and noisy) = %.6f\n', amax);
assert(amax*2^P.F_IN < 2^(P.W_IN-1), 'input exceeds the 16-bit range');
fprintf('    input: %d-bit two''s complement, %d fraction bits (S.%d), range [%g, %g), LSB = 2^-%d\n', ...
    P.W_IN, P.F_IN, P.F_IN, -2^(P.W_IN-1-P.F_IN), 2^(P.W_IN-1-P.F_IN), P.F_IN);
fprintf('    product Re/Im{z z*}: %d bits signed, |z|^2: %d bits unsigned\n', P.W_PROD, P.W_PWR);
fprintf('    accumulators (R = %d -> +%d bits): correlation %d bits signed, energy %d bits unsigned\n', ...
    R, log2(R), P.W_ACC_C, P.W_ACC_E);
fprintf('    normalizer keeps %d magnitude bits, numerator %d bits, denominator %d bits\n', P.K - 1, P.W_NUM, P.W_DEN);
fprintf('    divider: dividend %d bits (num << %d), divisor %d bits, quotient U2.%d (%d bits), eta = %d\n', ...
    P.W_NUM + P.F_PHI, P.F_PHI, P.W_DEN, P.F_PHI, P.W_PHI, P.ETA_Q);

% Bit-true model of the HDL and the test vectors of PATTERN.sv (02_hdl/00_TESTBED/src)
vecdir = fullfile(fileparts(mfilename('fullpath')), '..', '02_hdl', '00_TESTBED', 'src');
if ~exist(vecdir, 'dir'), mkdir(vecdir); end
fid = fopen(fullfile(vecdir, 'eta.hex'), 'w');
fprintf(fid, '%04X\n', P.ETA_Q);
fclose(fid);
o = cell(1, 2);
for c = 1:2
    o{c} = sbd_model(cases{c, 2}, P);
    err = o{c}.phi/2^P.F_PHI - min(cases{c, 3}, (2^P.W_PHI - 1)/2^P.F_PHI);
    fprintf('    [%s] max |Phi_q - Phi_float| = %.3g (= %.2f LSB of U2.14), max shift = %d, max |C| = 2^%.2f, max E = 2^%.2f\n', ...
        cases{c, 1}, max(abs(err)), max(abs(err))*2^P.F_PHI, max(o{c}.sh), ...
        log2(max(abs([o{c}.c_re; o{c}.c_im]))), log2(max(o{c}.e)));
    write_vectors(vecdir, cases{c, 1}, o{c}, P);
end

%% (d) Hardware results versus the floating-point simulation
% The VCS simulation (02_hdl/01_RTL/01_run) checks that the RTL output equals sbd_model bit by bit.
fprintf('\n(d) local peak searcher: eta = %.2f, decision %d samples after the candidate, flush at the last sample\n', P.ETA, P.W_WIN);
fh = new_fig('Q3(d)', [900 760]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    nmf = cases{c, 3};
    phq = o{c}.phi/2^P.F_PHI;
    det = o{c}.det;
    dphi = o{c}.det_phi/2^P.F_PHI;
    detf = sbd_peak_search(nmf, P.ETA, P.W_WIN);
    [~, kb] = arrayfun(@(s) max(nmf((s-1)*Ns + (1:Ns))), 1:nsym);
    kb = kb + (0:nsym-1)*Ns;
    fprintf('    [%s] HDL model : m = %s (Phi_q = %s)\n', cases{c, 1}, sprintf(' %d', det), sprintf(' %.4f', dphi));
    fprintf('    [%s] float, same searcher : m = %s\n', cases{c, 1}, sprintf(' %d', detf));
    fprintf('    [%s] Q2 per-160-sample-block argmax : m = %s\n', cases{c, 1}, sprintf(' %d', kb));

    nexttile;
    h1 = plot(1:M, nmf, ':', 'Color', p.muted, 'LineWidth', 1.4); hold on;
    h2 = plot(1:M, phq, 'Color', p.blue, 'LineWidth', 1.1);
    h3 = yline(P.ETA, '--', 'Color', p.red, 'LineWidth', 1.1, 'Layer', 'bottom');
    h4 = plot(kb, nmf(kb), 's', 'Color', p.orange, 'MarkerSize', 9, 'LineWidth', 1.4);
    xlim([1 M]); ylim([0 1.35]); xticks([1 40:40:M]);
    peak_labels(det, dphi, p.red, '%.4f');
    h5 = plot(nan, nan, 'o', 'Color', p.red, 'MarkerFaceColor', p.red, 'MarkerSize', 6);
    legend([h1 h2 h3 h4 h5], {'MATLAB float', 'HDL (U2.14)', sprintf('$\\eta = %.2f$', P.ETA), ...
        'Q2 per-block argmax', 'HDL detection'}, 'Location', 'north', 'Orientation', 'horizontal', 'FontSize', 10);
    ylabel('$\Phi_{NM}(m)$');
    title(sprintf('%s: HDL detections $m = %s$', cases{c, 4}, strjoin(compose('%d', det), ', ')));
end
tl_labels(tl, '', '$m$');
save_fig(fh, 'q3d');

%% Local functions
function v = lobe_peaks(phi, bnd, R, M)
% largest value of phi inside each CP correlation lobe |m - b| < R
v = arrayfun(@(b) max(phi(max(1, b - R + 1):min(M, b + R - 1))), bnd);
end

function write_vectors(vecdir, name, o, P)
% <name>_in.hex : one input sample per line, {Re, Im} as two 16-bit hex words
% <name>_phi.hex: golden divider output (U2.14) for every sample
% <name>_det.txt: number of detections, then "m phi_code" per detection
tc = @(x, w) mod(x, 2^w);                    % two's complement code
fid = fopen(fullfile(vecdir, [name '_in.hex']), 'w');
fprintf(fid, '%04X%04X\n', [tc(o.zr, P.W_IN) tc(o.zi, P.W_IN)].');
fclose(fid);
fid = fopen(fullfile(vecdir, [name '_phi.hex']), 'w');
fprintf(fid, '%04X\n', o.phi);
fclose(fid);
fid = fopen(fullfile(vecdir, [name '_det.txt']), 'w');
fprintf(fid, '%d\n', numel(o.det));
fprintf(fid, '%d %d\n', [o.det o.det_phi].');
fclose(fid);
end
