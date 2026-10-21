function res = polar_simulate(N, K, ebn0_list, L, crc_name, varargin)
%POLAR_SIMULATE  Monte-Carlo BER/BLER of a polar code, BPSK over AWGN.
%   K counts all bits on information positions (message + CRC).
%   Eb/N0 is per MESSAGE bit: rate = (K - r)/N.
%   Options (name/value): 'design', 2.0 | 'max_frames', 2e4 |
%                         'min_errors', 200 | 'batch', 500 | 'seed', 1
if nargin < 4 || isempty(L), L = 1; end
if nargin < 5, crc_name = ''; end
opt = struct('design', 2.0, 'max_frames', 2e4, 'min_errors', 200, ...
             'batch', 500, 'seed', 1);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end
rng(opt.seed);
crc = [];  r = 0;
if ~isempty(crc_name)
    r = crc_make(crc_name, 1).r;          % CRC length
    crc = crc_make(crc_name, K - r);
end
k = K - r;                                % message bits
rate = k / N;
info = ga_construct(N, K, opt.design, K/N);
res.ebn0 = ebn0_list;  res.ber = zeros(size(ebn0_list));
res.bler = res.ber;    res.errors = res.ber;
for s = 1:numel(ebn0_list)
    esn0 = 10^(ebn0_list(s)/10) * rate;
    sigma = sqrt(1/(2*esn0));
    be = 0; bl = 0; fr = 0;
    while fr < opt.max_frames && bl < opt.min_errors
        B = opt.batch;
        m = double(rand(B, k) < 0.5);
        if isempty(crc), pay = m; else, pay = crc_attach(crc, m); end
        u = zeros(B, N);  u(:, info) = pay;
        x = polar_encode(u);
        y = (1 - 2*x) + sigma * randn(B, N);          % BPSK + AWGN
        llr = 2*y / sigma^2;
        u_hat = polar_list_decode(llr, info, L, crc);
        m_hat = u_hat(:, info);  m_hat = m_hat(:, 1:k);
        e = (m_hat ~= m);
        be = be + sum(e(:));
        bl = bl + sum(any(e, 2));
        fr = fr + B;
    end
    res.ber(s) = be / (fr * k);
    res.bler(s) = bl / fr;
    res.errors(s) = bl;
    fprintf('  Eb/N0 = %4.1f dB   BER = %.2e   BLER = %.2e   (%d errors / %d frames)\n', ...
            ebn0_list(s), res.ber(s), res.bler(s), bl, fr);
end
end
