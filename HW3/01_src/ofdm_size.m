function [Ng, N] = ofdm_size(L)
% OFDM sizes of Q3 for an L-tap channel
% Ng: smallest power of 2 >= L - 1 (CP must cover the channel memory)
% N:  smallest power of 2 with Ng/N < 1/4, i.e. N > 4*Ng
Ng = 2^nextpow2(L - 1);
N = 2^nextpow2(4*Ng + 1);
end
