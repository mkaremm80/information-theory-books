%LAB5_LZW  LZW compressor and decompressor (MATLAB version of lzw.py).
c = lzw_encode('ABABABA', 'AB');
fprintf('ABABABA -> %s -> %s, %d bits\n', mat2str(c), lzw_decode(c, 'AB'), lzw_bits(c, 2));
t = 'TOBEORNOTTOBEORTOBEORNOT';
c = lzw_encode(t, '');
fprintf('%s\n  codes: %s\n  round trip OK: %d, %d bits (original %d)\n', t, mat2str(c), ...
        strcmp(lzw_decode(c, ''), t), lzw_bits(c, 256), 8*numel(t));
