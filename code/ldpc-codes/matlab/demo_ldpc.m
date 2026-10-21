% DEMO_LDPC  Worked examples of the book, computed with the MATLAB library.
H = [1 1 0 1 0 0; 0 1 1 0 1 0; 1 0 0 0 1 1; 0 0 1 1 0 1];
code = ldpc_from_H(H, 'toy');
fprintf('Toy code: rank %d, girth %d\n', gf2_rank(H), ldpc_girth(code));
[G, info] = generator_from_H(H);
fprintf('G = %s\n', mat2str(G)); fprintf('information positions (1-based): %s\n', mat2str(info));
fprintf('u = [1 0 1] -> c = %s\n', mat2str(encode_G(G, [1 0 1])));

% Chapter 7: peeling and bit flipping
y = [-1 0 -1 0 1 -1];
[x, rounds] = peel(code, y);
fprintf('peeling (?,0,?,0,1,?) -> %s in %d rounds\n', mat2str(x), rounds);
[x, it] = bit_flip(code, [0 0 1 1 1 1], 10);
fprintf('bit flipping 001111 -> %s in %d iteration(s)\n', mat2str(x), it);

% Chapter 8: sum-product example
L = 4 * [0.9 1.2 -0.8 -0.3 -1.1 -0.7];
c2v = check_update(code, L(code.vn), 'spa');
fprintf('SPA iteration-1 check messages: %s\n', mat2str(round(c2v * 1e4) / 1e4));
[x, it, app] = bp_decode(code, L, 20, 'spa');
fprintf('SPA: decoded %s after %d iteration(s); APP = %s\n', mat2str(x), it, mat2str(round(app * 1e4) / 1e4));

% Chapter 9: min-sum variants and layered decoding
for r = {'ms', 'nms', 'oms'}
    c2v = check_update(code, L(code.vn), r{1}, 0.75, 0.5);
    fprintf('%-3s check messages: %s\n', r{1}, mat2str(round(c2v * 1e4) / 1e4));
end
code = make_layers(code);
[x, it, app] = layered_decode(code, L, 20, 'spa');
fprintf('layered SPA: %s, %d iteration(s), APP = %s\n', mat2str(x), it, mat2str(round(app * 1e4) / 1e4));

% Chapter 10: thresholds
fprintf('BEC thresholds: (3,6) %.4f  (3,4) %.4f  (4,8) %.4f\n', threshold_bec([3 1], [6 1]), ...
        threshold_bec([3 1], [4 1]), threshold_bec([4 1], [8 1]));
fprintf('GA threshold (3,6): sigma = %.4f\n', threshold_ga([3 1], [6 1]));

% Chapter 11: 5G NR parameters
for KE = [120 360; 256 512; 1024 2048; 1024 1152; 8448 25344].'
    nr = nr_code(KE(1), KE(2));
    fprintf('K=%5d E=%5d: BG%d Z=%3d iLS=%d fillers=%3d rows used=%2d\n', KE(1), KE(2), nr.bg, nr.Z, ...
            nr.ils, nr.Kp - nr.K, nr.code.m / nr.Z);
end
