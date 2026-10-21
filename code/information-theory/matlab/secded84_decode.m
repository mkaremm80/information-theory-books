function [d, status] = secded84_decode(r)
%SECDED84_DECODE  Extended Hamming (8,4): r = [7-bit Hamming word, p0].
[~, ~, s] = hamming74_decode(r(1:7));
overall = mod(sum(r), 2);
if s == 0 && overall == 0
    d = r([3 5 6 7]); status = 'no error';
elseif overall == 1 && s == 0
    d = r([3 5 6 7]); status = 'parity bit error';
elseif overall == 1
    d = hamming74_decode(r(1:7)); status = 'corrected single error';
else
    d = []; status = 'double error detected';
end
end
