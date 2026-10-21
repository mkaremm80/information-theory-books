%FINAL_PROJECT  End-to-end polar-coded control link (MATLAB version of
%   final_project.py): A = 40 bits + CRC24C (K = 64), N = 128, BPSK/AWGN.
%   Receivers: SC (CRC used only for detection) and CA-SCL with L = 8.
%   Measures BER, BLER and undetected errors (wrong message passing the CRC).
N = 128;  A = 40;
crc  = crc_make('CRC24C', A);
K    = A + crc.r;                                 % 64 payload bits
R    = A / N;                                     % message rate 0.3125
info = ga_construct(N, K, 2.0, K/N);
names = {'SC', 'CA-SCL L=8'};  Ls = [1 8];  crcs = {[], crc};

% ---------------------------------------------------------- one frame
rng(2026);
m   = double(rand(1, A) < 0.5);                   % step 1: data
pay = crc_attach(crc, m);                         % step 2: CRC
u   = zeros(1, N);  u(info) = pay;                % step 3: map + encode
x   = polar_encode(u);
eb  = 4.5;  sigma = sqrt(1/(2*R*10^(eb/10)));
y   = (1 - 2*x) + sigma*randn(1, N);              % steps 4-5: BPSK + AWGN
llr = 2*y/sigma^2;
fprintf('hard-decision errors: %d of %d\n', sum((llr < 0) ~= (x == 1)), N);
for d = 1:2                                       % step 6: decode
    uh = polar_list_decode(llr, info, Ls(d), crcs{d});
    fprintf('%-11s: bit errors = %d, CRC pass = %d\n', names{d}, ...
            sum(uh(find(info, A)) ~= m), crc_check(crc, uh(info)));
end

% ---------------------------------------------------------- Monte-Carlo
E = 2:0.5:5.5;  B = 2000;
BLER = zeros(2, numel(E));  BER = BLER;  UND = BLER;
for d = 1:2
    fprintf('%s\n', names{d});
    rng(1);
    for s = 1:numel(E)
        sigma = sqrt(1/(2*R*10^(E(s)/10)));
        be = 0; bl = 0; und = 0; fr = 0;
        while fr < 1e5 && bl < 150
            m = double(rand(B, A) < 0.5);
            u = zeros(B, N);  u(:, info) = crc_attach(crc, m);
            llr = 2*((1 - 2*polar_encode(u)) + sigma*randn(B, N))/sigma^2;
            ph  = polar_list_decode(llr, info, Ls(d), crcs{d});
            ph  = ph(:, info);
            wrong = any(ph(:, 1:A) ~= m, 2);
            be  = be + sum(sum(ph(:, 1:A) ~= m));
            bl  = bl + sum(wrong);
            und = und + sum(wrong & crc_check(crc, ph));   % step 7
            fr  = fr + B;
        end
        BER(d, s) = be/(fr*A);  BLER(d, s) = bl/fr;  UND(d, s) = und;
        fprintf('  %4.1f dB  BER %.2e  BLER %.2e  undetected %d/%d\n', ...
                E(s), BER(d, s), BLER(d, s), und, fr);
    end
end
figure;
semilogy(E, BLER(1, :), 'o-', E, BLER(2, :), 's-', 'LineWidth', 1.6); grid on;
hold on; semilogy(E([1 end]), [1e-2 1e-2], 'k:');      % target BLER 10^-2
xlabel('E_b/N_0 (dB), message rate 40/128'); ylabel('BLER');
legend([names, {'target 10^{-2}'}], 'Location', 'southwest');
