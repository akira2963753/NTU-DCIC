function P = sbd_param()
% Fixed-point and algorithm parameters of the symbol boundary detector (Q3)
% Keep in sync with 01_src/02_RTL/define.vh
P.R = 32;                % correlation window
P.L = 128;               % correlation lag (= N, CP-based timing as in Q2(c))
P.W_IN = 16;             % input word length per I/Q (two's complement)
P.F_IN = 16;             % input fraction bits: z_q = round(z*2^F_IN), range [-0.5, 0.5)
P.W_PROD = 2*P.W_IN + 1;                 % z(m) z*(m-L) real/imag part: 33 bits signed
P.W_PWR = 2*P.W_IN;                      % |z(m)|^2: 32 bits unsigned
P.W_ACC_C = P.W_PROD + log2(P.R);        % correlation accumulator: 38 bits signed
P.W_ACC_E = P.W_PWR + log2(P.R);         % energy accumulator: 37 bits unsigned
P.K = 24;                % normalizer keeps K-1 = 23 magnitude bits before squaring
P.W_NUM = 2*(P.K - 1) + 1;               % |C|^2 after normalization: 47 bits
P.W_DEN = 2*(P.K - 1);                   % E^2 after normalization: 46 bits
P.W_PHI = 16;            % divider quotient, unsigned U2.14, saturates at 4 - 2^-14
P.F_PHI = 14;
P.ETA = 0.52;            % threshold of Q3(a): midpoint of the 99.9% quantiles (0.46 outside, 0.58 peak)
P.ETA_Q = round(P.ETA*2^P.F_PHI);
P.W_WIN = 32;            % peak searcher decides W_WIN samples after the candidate peak
end
