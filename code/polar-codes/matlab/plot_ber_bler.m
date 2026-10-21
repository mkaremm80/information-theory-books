%PLOT_BER_BLER  BER/BLER versus Eb/N0 for SC, SCL and CA-SCL (N=256, R=1/2).
%   Takes a few minutes in MATLAB; reduce 'max_frames' for a quick look.
N = 256;  E = 0.5:0.5:3.5;
opts = {'max_frames', 5e4, 'min_errors', 150, 'batch', 1000, 'seed', 2026};
disp('SC');            sc  = polar_simulate(N, 128, E, 1, '', opts{:});
disp('SCL L=8');       scl = polar_simulate(N, 128, E, 8, '', opts{:});
disp('CA-SCL L=8');    ca  = polar_simulate(N, 139, E, 8, 'CRC11', opts{:});
uncoded = 0.5 * erfc(sqrt(10.^(E/10)));                % Q(sqrt(2 Eb/N0))

figure('Position', [100 100 950 380]);
subplot(1, 2, 1);
semilogy(E, sc.bler, 'o-', E, scl.bler, 's-', E, ca.bler, '^-', 'LineWidth', 1.6);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BLER');
legend('SC', 'SCL, L=8', 'CA-SCL, L=8, CRC-11', 'Location', 'southwest');
title('N = 256, R = 1/2');
subplot(1, 2, 2);
semilogy(E, sc.ber, 'o-', E, scl.ber, 's-', E, ca.ber, '^-', E, uncoded, 'k:', 'LineWidth', 1.6);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('SC', 'SCL, L=8', 'CA-SCL, L=8', 'uncoded BPSK', 'Location', 'southwest');
title('N = 256, R = 1/2');
