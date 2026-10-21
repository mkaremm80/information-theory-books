%DEMO_POLAR  Small demonstrations with printed outputs (MATLAB version of
%            demo_ch11.py). Run from the folder that holds the library.
disp('1) Encoder -------------------------------------------');
info = logical([0 0 0 1 0 1 1 1]);             % A = {3,5,6,7} (0-based)
u = zeros(1, 8);  u(info) = [1 0 1 1];
x = polar_encode(u);
fprintf('u = %s   x = %s\n', mat2str(u), mat2str(x));
fprintf('matrix check: %d\n', isequal(x, mod(u * generator_matrix(8), 2)));

disp('2) Frozen-bit selection ------------------------------');
fprintf('BEC(0.5) Z, N=8 : %s\n', mat2str(round(1e4*bec_bhattacharyya(8, 0.5))/1e4));
[~, o] = sort(bec_bhattacharyya(8, 0.5));
fprintf('BEC info set K=4 (0-based): %s\n', mat2str(sort(o(1:4)) - 1));
fprintf('GA info set N=16, K=8 @2dB (0-based): %s\n', mat2str(find(ga_construct(16, 8, 2.0)) - 1));

disp('3) SC decoder -----------------------------------------');
llr1 = [-1.5 2.0 -0.8 1.2 -0.6 -2.5 1.8 -3.0];
fprintf('SC  u_hat = %s\n', mat2str(polar_list_decode(llr1, info, 1)));

disp('4) SCL decoder ----------------------------------------');
llr2 = [-2.0 0.3 -0.4 -2.0 1.0 -5.8 3.2 -4.8];
fprintf('SC  u_hat = %s (wrong)\n', mat2str(polar_list_decode(llr2, info, 1)));
[uh, ~, PM] = polar_list_decode(llr2, info, 2);
fprintf('SCL u_hat = %s   path metrics: %s\n', mat2str(uh), mat2str(PM));

disp('5) CRC ------------------------------------------------');
crc = crc_make('CRC6', 10);
fprintf('message + CRC6 : %s\n', mat2str(crc_attach(crc, [1 0 1 1 0 0 1 0 1 1])));

disp('6-7) BPSK + AWGN --------------------------------------');
rng(0);
esn0 = 10^(2.0/10) * 0.5;  sigma = sqrt(1/(2*esn0));
y = (1 - 2*x) + sigma * randn(1, 8);
fprintf('sigma = %.4f\n', sigma);
fprintf('y     = %s\n', mat2str(round(100*y)/100));
fprintf('LLR   = %s\n', mat2str(round(200*y/sigma^2)/100));

disp('8) BER / BLER -----------------------------------------');
res = polar_simulate(128, 64, [1 2 3], 1, '', 'max_frames', 2e4, ...
                     'min_errors', 100, 'batch', 1000, 'seed', 7);
