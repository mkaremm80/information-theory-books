%LAB7_HAMMING  Hamming (7,4) and SECDED (MATLAB version of hamming.py).
d = [1 0 1 1];
c = hamming74_encode(d);
fprintf('data %s -> codeword %s\n', mat2str(d), mat2str(c));
r = c; r(5) = 1 - r(5);                          % error at position 5
[dd, fixed, s] = hamming74_decode(r);
fprintf('received %s syndrome %d -> corrected %s data %s\n', mat2str(r), s, mat2str(fixed), mat2str(dd));
r2 = c; r2([2 6]) = 1 - r2([2 6]);               % two errors
[dd2, ~, s2] = hamming74_decode(r2);
fprintf('two errors -> syndrome %d data %s (WRONG)\n', s2, mat2str(dd2));
c8 = [c mod(sum(c), 2)];  r8 = c8; r8([2 6]) = 1 - r8([2 6]);
[~, st] = secded84_decode(r8);
fprintf('SECDED with the same two errors -> %s\n', st);
