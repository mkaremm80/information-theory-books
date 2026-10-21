%TEST_POLAR_MATLAB  Self-tests of the MATLAB polar library.
% 1) encoder equals u*G_N
for N = [2 4 8 16 64 256]
    u = double(rand(50, N) < 0.5);
    assert(isequal(polar_encode(u), mod(u * generator_matrix(N), 2)));
end
disp('encoder ok');
% 2) worked examples of the book
info = logical([0 0 0 1 0 1 1 1]);
assert(isequal(polar_encode([0 0 0 1 0 0 1 1]), [1 0 1 0 0 1 0 1]));
assert(isequal(polar_list_decode([-1.5 2.0 -0.8 1.2 -0.6 -2.5 1.8 -3.0], info, 1), [0 0 0 1 0 0 1 1]));
llr2 = [-2.0 0.3 -0.4 -2.0 1.0 -5.8 3.2 -4.8];
assert(isequal(polar_list_decode(llr2, info, 1), [0 0 0 0 0 0 1 1]));
[uh, ~, PM] = polar_list_decode(llr2, info, 2);
assert(isequal(uh, [0 0 0 1 0 0 1 1]) && max(abs(PM - [2.0 2.7])) < 1e-12);
disp('book examples ok');
% 3) CRC
assert(isequal(crc_bits([1 1 0 1], [3 1 0]), [0 0 1]));
c = crc_make('CRC6', 10);
w = crc_attach(c, [1 0 1 1 0 0 1 0 1 1]);
assert(isequal(w(11:16), [0 1 1 0 1 1]) && crc_check(c, w));
w(3) = 1 - w(3);  assert(~crc_check(c, w));
disp('crc ok');
% 4) identical decisions to the Python library on random noisy blocks
%    (pyref.mat holds LLRs and decisions produced by the Python code)
R = load('pyref.mat');
tags = {'N64L1', 'N64L4', 'N128L8c', 'N256L8'};
Ls = [1 4 8 8];
for t = 1:numel(tags)
    llr = R.([tags{t} '_llr']);  info = logical(R.([tags{t} '_info']));
    if tags{t}(end) == 'c'
        crc = crc_make('CRC11', sum(info) - 11);
    else
        crc = [];
    end
    uh = polar_list_decode(llr, info, Ls(t), crc);
    fprintf('%-8s same decisions as Python: %d of %d blocks\n', tags{t}, ...
            sum(all(uh == R.([tags{t} '_uh']), 2)), size(llr, 1));
end
% 5) GA construction equals Python
[info, m] = ga_construct(1024, 512, 2.0, 0.5);
p = R.ga1024;  big = p > 1e-3;           % ignore near-zero means (useless channels)
[~, o] = sort(p);  info_py = false(1, 1024);  info_py(o(end-511:end)) = true;
fprintf('GA: max relative difference %.1e, same information set: %d\n', ...
        max(abs(m(big) - p(big)) ./ p(big)), isequal(info, info_py));
