function d = hw4_data()
% Data of Q2 and Q3: OFDMTx from HW4.mat, frame parameters and the 12 dB noisy record
s = load(fullfile(fileparts(mfilename('fullpath')), '..', 'HW4.mat'));
d.x = s.OFDMTx(:);       % 640 x 1: preamble + 3 normal OFDM symbols
d.N = 128;               % FFT size
d.Ng = 32;               % cyclic prefix length
d.Ns = d.N + d.Ng;       % samples per OFDM symbol
d.nsym = numel(d.x)/d.Ns;

% Record of Q2(c)/(d) and Q3: preamble discarded, 480 samples
d.z = d.x(d.Ns+1:end);

% Circular complex white Gaussian noise, SNR = mean(|z|^2)/sigma2 over the 480 samples
d.snr_db = 12;
d.seed = 11407439;       % the same noisy record is used in Q2(d) and Q3
d.sigma2 = mean(abs(d.z).^2)/10^(d.snr_db/10);
d.noise = @() sqrt(d.sigma2/2)*(randn(size(d.z)) + 1j*randn(size(d.z)));
rng(d.seed);
d.zn = d.z + d.noise();
end
