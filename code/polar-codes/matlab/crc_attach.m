function w = crc_attach(crc, m)
%CRC_ATTACH  Append CRC bits to every row of m:  w = [m | mod(m*P,2)].
w = [m, mod(m * crc.P, 2)];
end
