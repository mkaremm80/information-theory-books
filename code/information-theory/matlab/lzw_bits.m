function total = lzw_bits(codes, n0)
%LZW_BITS  Bits with dynamic widths: code k (0-based) uses ceil(log2(n0+k)).
k = 0:numel(codes) - 1;
total = sum(max(1, ceil(log2(n0 + k))));
end
