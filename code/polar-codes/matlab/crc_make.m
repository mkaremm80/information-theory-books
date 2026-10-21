function crc = crc_make(name, k)
%CRC_MAKE  Build a CRC as a linear map: crc bits = mod(m * P, 2).
%   name : 'CRC6', 'CRC11' or 'CRC24C' (3GPP TS 38.212, clause 5.1)
%   k    : message length in bits.
switch upper(name)
    case 'CRC6',   e = [6 5 0];
    case 'CRC11',  e = [11 10 9 5 0];
    case 'CRC24C', e = [24 23 21 20 17 15 13 12 8 4 2 1 0];
    otherwise, error('unknown CRC');
end
crc.r = e(1);
crc.k = k;
crc.exps = e;
crc.P = zeros(k, crc.r);
for i = 1:k
    unit = zeros(1, k);  unit(i) = 1;
    crc.P(i, :) = crc_bits(unit, e);
end
end
