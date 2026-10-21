%COMPARE_LENGTHS_RATES  SC decoding: effect of block length and code rate.
opts = {'max_frames', 5e4, 'min_errors', 150, 'batch', 1000, 'seed', 2026};
figure('Position', [100 100 950 380]);
subplot(1, 2, 1); hold on;
E = 0:0.5:4.5;
for N = [64 128 256 512 1024]
    fprintf('N = %d\n', N);
    r = polar_simulate(N, N/2, E, 1, '', opts{:});
    ok = r.errors >= 20;                       % plot reliable points only
    semilogy(E(ok), r.bler(ok), 'o-', 'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', N));
end
set(gca, 'YScale', 'log'); grid on; legend('Location', 'southwest');
xlabel('E_b/N_0 (dB)'); ylabel('BLER'); title('SC, R = 1/2');
subplot(1, 2, 2); hold on;
E = 0:0.5:5.5;
for K = [64 128 192]
    fprintf('K = %d\n', K);
    r = polar_simulate(256, K, E, 1, '', opts{:});
    ok = r.errors >= 20;
    semilogy(E(ok), r.bler(ok), 'o-', 'LineWidth', 1.5, 'DisplayName', sprintf('R = %d/256', K));
end
set(gca, 'YScale', 'log'); grid on; legend('Location', 'southwest');
xlabel('E_b/N_0 (dB)'); ylabel('BLER'); title('SC, N = 256');
