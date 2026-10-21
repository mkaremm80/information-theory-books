function ok = crc_check(crc, w)
%CRC_CHECK  True for every row of w = [m | p] whose CRC is correct.
m = w(:, 1:crc.k);
p = w(:, crc.k+1:end);
ok = all(mod(m * crc.P, 2) == p, 2);
end
