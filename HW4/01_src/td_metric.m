function [phi_dc, phi_nm, num, den] = td_metric(z, R, L)
% Delay-and-correlate and normalized timing metrics for m = 1..numel(z)
%   phi_dc(m) = | sum_{r=0}^{R-1} z(m-r) z*(m-r-L) |
%   phi_nm(m) = phi_dc(m)^2 / ( sum_{r=0}^{R-1} |z(m-r)|^2 )^2
% with z(m) = 0 for m <= 0 (zeros are prepended, nothing is wrapped)
z = z(:);
M = numel(z);
zp = [zeros(R + L, 1); z];          % zp(k + R + L) = z(k)
off = R + L;
c = zeros(M, 1);
e = zeros(M, 1);
for m = 1:M
    idx = off + m - (0:R-1);        % m - r
    c(m) = sum(zp(idx).*conj(zp(idx - L)));
    e(m) = sum(abs(zp(idx)).^2);
end
phi_dc = abs(c);
num = abs(c).^2;
den = e.^2;
phi_nm = num./den;
end
