% SIM_LDPC_MATLAB  FER of the 5G NR code K = 1024, R = 1/2 (BG2, Z = 104) with
% layered normalised min-sum (alpha = 0.75, 50 iterations) in MATLAB/Octave.
rng_seed = 2026; randn('seed', rng_seed); rand('seed', rng_seed); %#ok<RAND>
nr = nr_code(1024, 2048);
nr.code = make_layers(nr.code, nr.Z);
enc = @(u) nr_encode(nr, u);
dec = @(L) layered_decode(nr.code, nr_llr_in(nr, L), 50, 'nms');
res = ldpc_simulate(enc, dec, 1024, 0.5, 0.5:0.25:2.0, 1:1024, 4000, 50, 200);
save('-v7', 'sim_matlab.mat', 'res');
