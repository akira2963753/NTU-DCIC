%% HW4 Q1: ULA channel and DFT receive beamforming
clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'plot_util'));

%% Parameters
Nr = 64;                 % receive antennas (= N of the DFT codebook)
N = Nr;
L = 2;                   % number of paths
kd = pi;                 % k*d (half-wavelength spacing)
theta_t1 = pi/7;         % AoD of path 1 (single Tx antenna: no Tx steering vector)
sid = 'M11407439';
m1 = sid(end) - '0';     % last digit of the student ID
m = -N/2:N/2-1;          % codebook indices

p = hw_palette();
rng(11407439);           % fixed seed so that the reported values are reproducible

a_ula = @(th) exp(-1j*kd*(0:Nr-1).'*sin(th(:).'))/sqrt(Nr);    % Nr x numel(th)
c_dft = @(mm) exp(-1j*2*pi*(0:N-1).'*mm(:).'/N)/sqrt(Nr);       % Nr x numel(mm)
phi_nom = @(mm) asin(2*mm/N);                                   % nominal pointing angle of c(m)

%% (a) AoA, path gains and channel vector
theta_r = (2*rand(1, L) - 1)*pi/3;         % uniform in (-pi/3, pi/3)
amag = [0.5 + 0.5*rand, 0.5*rand];         % |alpha_1| > 0.5 > |alpha_2|
aph = (2*rand(1, L) - 1)*pi;               % uniform phase in (-pi, pi)
alpha = amag.*exp(1j*aph);
h = a_ula(theta_r)*alpha.'/sqrt(L);        % 64 x 1

assert(all(abs(theta_r) < pi/3) && abs(alpha(1)) > 0.5 && abs(alpha(2)) < 0.5);
fprintf('(a) Nr = %d, kd = pi, L = %d, theta_t,1 = pi/7 = %.6f rad = %.4f deg\n', Nr, L, theta_t1, rad2deg(theta_t1));
for l = 1:L
    fprintf('    theta_r,%d = %+.6f rad = %+.4f deg\n', l, theta_r(l), rad2deg(theta_r(l)));
end
for l = 1:L
    fprintf('    alpha_%d = %+.6f %+.6fj, |alpha_%d| = %.6f, angle = %+.6f rad = %+.4f deg\n', ...
        l, real(alpha(l)), imag(alpha(l)), l, amag(l), aph(l), rad2deg(aph(l)));
end
fprintf('    size(h) = %d x %d, ||h||^2 = %.6f\n', size(h, 1), size(h, 2), norm(h)^2);

%% (b) Beam pattern of c(m1)
th = linspace(-pi/2, pi/2, 4001);
bp = abs(c_dft(m1)'*a_ula(th));
[bpmax, ibp] = max(bp);
fprintf('\n(b) m1 = %d (last digit of %s)\n', m1, sid);
fprintf('    nominal pointing angle phi(m1) = asin(2*m1/N) = %+.6f rad = %+.4f deg\n', phi_nom(m1), rad2deg(phi_nom(m1)));
fprintf('    scanned peak |c(m1)^H a(theta)| = %.6f at theta = %+.4f deg\n', bpmax, rad2deg(th(ibp)));
fprintf('    first nulls at sin(theta) = 2*m1/N +- 2/N -> theta = %+.4f, %+.4f deg\n', ...
    rad2deg(asin(2*(m1 - 1)/N)), rad2deg(asin(2*(m1 + 1)/N)));

fh = new_fig('Q1(b)', [900 600]);
tl = tiledlayout(1, 1, 'Padding', 'compact');
pax = polaraxes(tl);
polarplot(pax, th, bp, 'Color', p.blue, 'LineWidth', 1.6); hold(pax, 'on');
polarplot(pax, [0 phi_nom(m1)], [0 1], '--', 'Color', p.red, 'LineWidth', 1.2);
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaLim = [-90 90];
pax.RLim = [0 1];
pax.ThetaTick = -90:15:90;
pax.TickLabelInterpreter = 'latex';
pax.FontSize = 12;
pax.GridColor = p.grid; pax.GridAlpha = 1;
legend(pax, {'$|\mathbf{c}(m_1)^H\mathbf{a}^{ULA}(\theta)|$', ...
    sprintf('$\\varphi(m_1) = %.2f^\\circ$', rad2deg(phi_nom(m1)))}, 'Location', 'southoutside');
tl_labels(tl, sprintf('Beam pattern of $\\mathbf{c}(m_1)$, $m_1 = %d$, $N_r = %d$', m1, Nr), '');
save_fig(fh, 'q1b');

%% (c) Combiner output versus m
g = abs(c_dft(m)'*h).';
fh = new_fig('Q1(c)', [900 460]);
stem(m, g, 'filled', 'Color', p.blue, 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
[gmax, imax] = max(g);
stem(m(imax), gmax, 'filled', 'Color', p.red, 'MarkerSize', 7, 'LineWidth', 2);
text(m(imax) + 1, gmax, sprintf('$\\hat m = %d$, $|\\mathbf{c}(\\hat m)^H\\mathbf{h}| = %.4f$', m(imax), gmax), ...
    'Color', p.red, 'FontSize', 12);
xlim([-33 32]); xticks(-32:8:32);
set(gca, 'XTickLabelRotation', 0, 'XMinorGrid', 'on');
ax = gca; ax.XAxis.MinorTickValues = -32:2:32;
xlabel('Beam index $m$');
ylabel('$|\mathbf{c}(m)^H\mathbf{h}|$');
title('Magnitude of the combiner output');
save_fig(fh, 'q1c');

fprintf('\n(c) |c(m)^H h| for m = -32 ~ 31:\n');
for k = 1:8:N
    fprintf('    m = %+3d ~ %+3d: %s\n', m(k), m(k + 7), sprintf(' %.4f', g(k:k+7)));
end

%% (d) Selected beam and pointing error
mh = m(imax);
phi = phi_nom(mh);
fprintf('\n(d) m_hat = argmax |c(m)^H h| = %d, |c(m_hat)^H h| = %.6f\n', mh, gmax);
fprintf('    phi(m_hat) = asin(2*m_hat/N) = %+.6f rad = %+.4f deg\n', phi, rad2deg(phi));
for l = 1:L
    fprintf('    |phi(m_hat) - theta_r,%d| = %.6f rad = %.4f deg\n', l, abs(phi - theta_r(l)), rad2deg(abs(phi - theta_r(l))));
end

%% (e) Supporting values
u = N*sin(theta_r)/2;                       % path direction in units of beam index
comp = (alpha/sqrt(L)).*(c_dft(mh)'*a_ula(theta_r));    % per-path contribution at m_hat
fprintf('\n(e) angular separation |theta_r,1 - theta_r,2| = %.4f deg, beam spacing in sin(theta) = 2/N = %.4f\n', ...
    rad2deg(abs(diff(theta_r))), 2/N);
for l = 1:L
    ml = round(u(l));
    fprintf(['    path %d: |alpha| = %.4f, angle(alpha) = %+.4f deg, N*sin(theta)/2 = %+.4f -> nearest beam %+d ' ...
        '(|c^H h| = %.4f); at m_hat: |c^H a(theta)| = %.4f, contribution = %.4f, angle = %+.4f deg\n'], ...
        l, amag(l), rad2deg(aph(l)), u(l), ml, g(m == ml), abs(c_dft(mh)'*a_ula(theta_r(l))), abs(comp(l)), rad2deg(angle(comp(l))));
end
fprintf('    |sum of contributions| = %.6f (check |c(m_hat)^H h| = %.6f)\n', abs(sum(comp)), gmax);
