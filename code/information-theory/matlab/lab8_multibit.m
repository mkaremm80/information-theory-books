%LAB8_MULTIBIT  (5,1) repetition and BCH(15,7) (MATLAB version of Lab 8).
msg = [1 0 1];
c = repelem(msg, 5);
r = c; r([2 4 7 10]) = 1 - r([2 4 7 10]);       % 2 errors in blocks 1 and 2
[dec, votes] = rep_decode(r);
fprintf('repetition: votes %s -> decoded %s\n', mat2str(votes), mat2str(dec));
r3 = c; r3(1:3) = 1 - r3(1:3);
fprintf('3 errors in one block -> %s (failure expected)\n', mat2str(rep_decode(r3)));
m = [1 0 1 1 0 0 1];
cw = bch1517_encode(m);
fprintf('BCH codeword %s\n', mat2str(cw));
rb = cw; rb([3 12]) = 1 - rb([3 12]);           % errors at 0-based positions 2 and 11
[mh, ~, pos] = bch1517_decode(rb);
fprintf('errors found at positions %s, message OK: %d\n', mat2str(pos), isequal(mh, m));
