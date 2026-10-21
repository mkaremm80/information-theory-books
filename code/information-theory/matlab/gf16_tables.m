function [EXP, LOG] = gf16_tables()
%GF16_TABLES  EXP(i+1) = alpha^i (i = 0..29) as integers 0..15, from x^4+x+1.
%   LOG(v+1) = i with alpha^i = v (v = 1..15).
EXP = zeros(1, 30); LOG = zeros(1, 16); v = 1;
for i = 0:14
    EXP(i+1) = v; LOG(v+1) = i;
    v = v * 2;
    if v >= 16, v = bitxor(v, 19); end          % 19 = 10011b
end
EXP(16:30) = EXP(1:15);
end
