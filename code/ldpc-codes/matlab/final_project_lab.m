% FINAL_PROJECT_LAB  MATLAB version of the final project (Chapter 16): one 5G NR
% transport block with TB CRC, segmentation, LDPC coding, two redundancy
% versions (HARQ-IR) and layered decoding.
A = 1000; R = 1/2; ebn0 = 0.75:0.25:1.75; frames = 600; batch = 100;
L_tb = 16 + 8 * (A > 3824);                          % CRC16 or CRC24A
if A > 3824, poly = [24 23 18 17 14 11 10 7 6 5 4 3 1 0]; else, poly = [16 12 5 0]; end
P = crc_matrix(poly, A);
B = A + L_tb;
if A <= 292 || (A <= 3824 && R <= 0.67) || R <= 0.25, bg = 2; else, bg = 1; end
nr = nr_code(B, 1e6, bg);                            % mother code (single code block here)
nr.code = make_layers(nr.code, nr.Z);
Z = nr.Z;
if bg == 1, N = 66 * Z; f = [0 17 33 56]; else, N = 50 * Z; f = [0 13 25 43]; end
buf = 2 * Z + 1:2 * Z + N; buf = buf(~ismember(buf, nr.filler));
G = round(A / R);
pos = @(rv) buf(mod(find(buf >= 2 * Z + 1 + f(rv + 1) * Z, 1) - 1 + (0:G - 1), numel(buf)) + 1);
fprintf('A = %d, R = %.2f: TB CRC %d bits, BG%d, Z = %d, fillers %d, k0 = %s\n', A, R, L_tb, bg, Z, ...
        nr.Kp - B, mat2str(f * Z));
randn('seed', 7); rand('seed', 7); %#ok<RAND>
for eb = ebn0
    e1 = 0; e2 = 0; n = 0;
    while n < frames
        a = randi([0 1], batch, A);
        b = [a, mod(a * P, 2)];                      % attach the TB CRC
        c = nr_encode(nr, b);
        soft = zeros(batch, nr.n);
        for t = 1:2                                  % rv0, then rv2
            p = pos(2 * (t - 1));
            soft(:, p) = soft(:, p) + awgn_llr(c(:, p), eb, R);
            Ld = soft; Ld(:, nr.filler) = 1e3;
            x = layered_decode(nr.code, Ld, 25, 'nms');
            bh = x(:, 1:B);
            ok = all(mod(bh(:, 1:A) * P, 2) == bh(:, A + 1:B), 2);
            if t == 1, e1 = e1 + sum(~ok); else, e2 = e2 + sum(~ok); end
        end
        n = n + batch;
    end
    fprintf('Eb/N0 = %4.2f dB: TB error after rv0 %.4f, after rv0+rv2 %.4f\n', eb, e1 / n, e2 / n);
end
