function ch = hw3_channel()
% Four-path channel shared by Q1 ~ Q4: PDP of Q1 and the seed-39 path gains of Q2(a)
m = 39;    % last two digits of the student ID

ch.m = m;
ch.P_dB = [0 -2 -6 -9];
ch.tau = [0 16 92 200]*1e-9;
ch.Fs = 250e6;

% Q1(a): linear and normalized powers
ch.p = 10.^(ch.P_dB/10);
ch.pt = ch.p/sum(ch.p);

% Q1(d): integer sample delays n_i = tau_i*Fs and discrete channel length
ch.ni = round(ch.tau*ch.Fs);
ch.L = max(ch.ni) + 1;

% Q2(a): gamma_i ~ CN(0,1), real and imaginary parts each with variance 1/2
rng(m, 'twister');
ch.gamma = (randn(1, 4) + 1j*randn(1, 4))/sqrt(2);
ch.g = sqrt(ch.p).*ch.gamma;
ch.h = ch.g/sqrt(sum(abs(ch.g).^2));

% Q4(b): h[n] = sum_i h_i delta[n - n_i], n = 0 ~ L-1
ch.hn = zeros(1, ch.L);
ch.hn(ch.ni + 1) = ch.h;
end
