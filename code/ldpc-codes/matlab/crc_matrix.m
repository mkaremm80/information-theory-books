function P = crc_matrix(poly, k)
%CRC_MATRIX crc(m) = mod(m*P, 2) for a k-bit message. poly lists the exponents
% of the generator, e.g. [16 12 5 0] for CRC16 of TS 38.212.
r = poly(1);
g = zeros(1, r + 1); g(r - poly + 1) = 1;            % g(1) = coefficient of x^r
P = zeros(k, r);
cur = g(2:end);                                      % x^r mod g
for i = k:-1:1
    P(i, :) = cur;
    top = cur(1);
    cur = [cur(2:end), 0];                           % multiply by x
    if top, cur = mod(cur + g(2:end), 2); end
end
end
