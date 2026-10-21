% TEST_LDPC_MATLAB  Compare the MATLAB library with the Python reference results.
ref = load('pyref_ldpc.mat');
d = @(x) double(x);
npass = 0;
report = @(name, ok) fprintf('%-62s %s\n', name, pass_fail(ok));

% 1. GF(2) and generator matrix
H6 = [1 1 0 1 0 0; 0 1 1 0 1 0; 1 0 0 0 1 1; 0 0 1 1 0 1];
[G6, info6] = generator_from_H(H6);
ok = gf2_rank(H6) == 3 && isequal(G6, [0 1 1 1 0 0; 1 1 0 0 1 0; 1 1 1 0 0 1]) && isequal(info6, [4 5 6]);
report('GF(2): toy code rank, G and information positions', ok); npass = npass + ok;
H = d(ref.H_peg); code = ldpc_from_H(H, 'PEG');
[G, info] = generator_from_H(H);
ok = isequal(G, d(ref.G_peg)) && isequal(info - 1, d(ref.info_peg)) && ~any(any(mod(H * G.', 2)));
report('GF(2): generator matrix of PEG(1008,504) identical to Python', ok); npass = npass + ok;

% 2. structure
ok = ldpc_girth(code) == 8 && count_4cycles(code) == 0 && ldpc_girth(ldpc_from_H(H6)) == 6;
q = qc_lift([0 0 0 0; 0 1 1 3], 5);
ok = ok && ldpc_girth(q) == 4 && count_4cycles(q) == 5 && ldpc_girth(qc_lift([0 0 0 0; 0 1 2 3], 5)) == 8;
report('Girth and 4-cycles (PEG 8, toy 6, QC 4 and 8)', ok); npass = npass + ok;

% 3. check-node update on random inputs
ok = true;
rules = {'spa', 'ms', 'nms', 'oms'};
for r = 1:4
    out = check_update(code, ref.v2c, rules{r}, 0.75, 0.5);
    ok = ok && max(abs(out(:) - ref.(['cu_' rules{r}])(:))) < 1e-9;
end
report('Check-node update (SPA, MS, NMS, OMS) equals Python', ok); npass = npass + ok;

% 4. flooding decoders on 200 noisy frames
L = ref.L_peg; c = d(ref.c_peg); ok = true; agree = 0;
for r = 1:4
    [x, it] = bp_decode(code, L, 50, rules{r});
    same = all(x == d(ref.(['x_' rules{r}])), 2) & it == reshape(d(ref.(['it_' rules{r}])), [], 1);
    agree = agree + sum(same);
end
report(sprintf('Flooding decoders: %d of 800 frames identical (decisions and iterations)', agree), agree == 800);
npass = npass + (agree == 800);

% 5. layered decoders
code = make_layers(code); agree = 0;
for r = {'nms', 'spa'}
    [x, it] = layered_decode(code, L, 50, r{1});
    agree = agree + sum(all(x == d(ref.(['xl_' r{1}])), 2) & it == reshape(d(ref.(['itl_' r{1}])), [], 1));
end
report(sprintf('Layered decoders (8 layers): %d of 400 frames identical', agree), agree == 400);
npass = npass + (agree == 400);

% 6. bit flipping and peeling
[x, it] = bit_flip(code, d(ref.r_bf), 50);
ok = isequal(x, d(ref.x_bf)) && isequal(it, reshape(d(ref.it_bf), [], 1));
[y, rd] = peel(code, d(ref.y_peel));
ok = ok && isequal(y, d(ref.x_peel)) && isequal(rd, reshape(d(ref.rd_peel), [], 1));
report('Bit flipping (50 frames) and peeling (50 frames) identical', ok); npass = npass + ok;

% 7. 5G NR codes: parameters and codewords
ok = true;
for i = 0:7
    KE = d(ref.(sprintf('nr%d_KE', i)));
    nr = nr_code(KE(1), KE(2));
    [cc, xx] = nr_encode(nr, d(ref.(sprintf('nr%d_u', i))));
    par = [nr.bg nr.Z nr.ils nr.code.m nr.code.n nr.code.E];
    ok = ok && isequal(par, d(ref.(sprintf('nr%d_par', i)))) && isequal(cc, d(ref.(sprintf('nr%d_c', i)))) ...
            && isequal(xx, d(ref.(sprintf('nr%d_x', i)))) && ~any(any(ldpc_syndrome(nr.code, cc)));
end
report('5G NR: parameters and codewords identical for 8 (K,E) pairs', ok); npass = npass + ok;

% 8. 5G NR layered decoding
nr = nr_code(1024, 2048); nr.code = make_layers(nr.code, nr.Z);
[xh, it] = layered_decode(nr.code, nr_llr_in(nr, ref.nr_L), 50, 'nms');
agree = sum(all(xh == d(ref.nr_xh), 2) & it == reshape(d(ref.nr_it), [], 1));
report(sprintf('5G NR layered NMS (K=1024, R=1/2): %d of 100 frames identical', agree), agree == 100);
npass = npass + (agree == 100);

% 9. density evolution
t1 = threshold_bec([3 1], [6 1]); t2 = threshold_bec([3 1], [4 1]); t3 = threshold_ga([3 1], [6 1]);
ok = abs(t1 - ref.th_bec36) < 1e-9 && abs(t2 - ref.th_bec34) < 1e-9 && abs(t3 - ref.th_ga36) < 2e-4;
report(sprintf('Thresholds: BEC %.4f %.4f, GA sigma %.4f', t1, t2, t3), ok); npass = npass + ok;
% 11. constructions that use MATLAB's own random stream
rand('seed', 3); randn('seed', 3); %#ok<RAND>
pc = peg_code(96, 48, 3); Hg = gallager_H(60, 3, 6);
ok = ldpc_girth(pc) >= 6 && all(pc.dv == 3) && gf2_rank(Hg) <= 30 - 2 && all(sum(Hg, 1) == 3) && all(sum(Hg, 2) == 6);
report(sprintf('PEG(96,48): girth %d; Gallager (3,6), n=60: rank %d', ldpc_girth(pc), gf2_rank(Hg)), ok);
npass = npass + ok;
fprintf('\n%d of 11 tests passed\n', npass);

