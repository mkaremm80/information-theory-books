function p = crc_bits(msg, exps)
%CRC_BITS  CRC remainder of msg(D)*D^r divided by g(D), by long division.
%   exps: exponents of g(D) with coefficient 1, highest first, e.g. [3 1 0].
r = exps(1);
g = zeros(1, r+1);
g(r - exps + 1) = 1;                       % g(1) is the D^r coefficient
reg = [msg(:).' zeros(1, r)];
for i = 1:numel(msg)
    if reg(i)
        reg(i:i+r) = xor(reg(i:i+r), g);
    end
end
p = reg(end-r+1:end);
end
