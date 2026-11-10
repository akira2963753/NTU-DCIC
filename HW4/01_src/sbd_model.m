function o = sbd_model(z, P)
% Bit-true model of the HDL symbol boundary detector (02_hdl/01_RTL)
%   z      : complex input record, m = 1..M (quantized here to P.W_IN bits per I/Q)
%   o.zr, o.zi : quantized input codes (two's complement integers), the HDL input
%   o.phi  : divider output Phi_NM in U2.14 integer code, one per input sample
%   o.det  : detected boundaries (sample index m), o.det_phi their metric codes
% All intermediate values stay below 2^53 (doubles are exact) except the divider
% dividend, which is handled in uint64.
qmax = 2^(P.W_IN-1);
quant = @(x) min(max(round(x(:)*2^P.F_IN), -qmax), qmax - 1);
zr = quant(real(z)); zi = quant(imag(z));
o.zr = zr; o.zi = zi;
M = numel(zr);
R = P.R; L = P.L;

% Delay correlator and energy accumulator (recursive form, z = 0 before m = 1)
zrL = [zeros(L, 1); zr(1:end-L)];
ziL = [zeros(L, 1); zi(1:end-L)];
pr = zr.*zrL + zi.*ziL;              % Re{z(m) z*(m-L)}
pim = zi.*zrL - zr.*ziL;             % Im{z(m) z*(m-L)}
pe = zr.^2 + zi.^2;                  % |z(m)|^2
prR = [zeros(R, 1); pr(1:end-R)];
piR = [zeros(R, 1); pim(1:end-R)];
peR = [zeros(R, 1); pe(1:end-R)];
o.c_re = cumsum(pr - prR);
o.c_im = cumsum(pim - piR);
o.e = cumsum(pe - peR);
assert(all(abs(pr) < 2^(P.W_PROD-1)) && all(abs(o.c_re) < 2^(P.W_ACC_C-1)) && all(abs(o.c_im) < 2^(P.W_ACC_C-1)));
assert(all(o.e < 2^P.W_ACC_E) && all(o.e >= 0));

% Block normalization: common right shift so that every magnitude fits in K-1 bits
mre = abs(o.c_re); mim = abs(o.c_im); me = o.e;
v = bitor(bitor(uint64(mre), uint64(mim)), uint64(me));
msb = zeros(M, 1);
nz = v > 0;
msb(nz) = floor(log2(double(v(nz))));          % exact for integers < 2^53
o.sh = max(0, msb - (P.K - 2));
o.sh(~nz) = 0;
mre = floor(mre./2.^o.sh); mim = floor(mim./2.^o.sh); me = floor(me./2.^o.sh);
assert(all([mre; mim; me] < 2^(P.K-1)));

% Squares and divider: phi = min(floor(num*2^14/den), 2^16-1), phi = 0 when den = 0
o.num = mre.^2 + mim.^2;
o.den = me.^2;
assert(all(o.num < 2^P.W_NUM) && all(o.den < 2^P.W_DEN));
q = zeros(M, 1);
k = o.den > 0;
q(k) = double(idivide(uint64(o.num(k))*uint64(2^P.F_PHI), uint64(o.den(k)), 'floor'));
o.phi = min(q, 2^P.W_PHI - 1);

% Threshold judgement and local peak searcher
[o.det, o.det_phi] = sbd_peak_search(o.phi, P.ETA_Q, P.W_WIN);
end
