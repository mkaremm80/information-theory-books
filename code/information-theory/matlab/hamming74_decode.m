function [d, r, s] = hamming74_decode(r)
%HAMMING74_DECODE  Single-error correction. s = syndrome = error position.
H = [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1];   % rows s1, s2, s4
sv = mod(H * r(:), 2);
s = sv(1) + 2*sv(2) + 4*sv(3);
if s > 0, r(s) = 1 - r(s); end
d = r([3 5 6 7]);
end
