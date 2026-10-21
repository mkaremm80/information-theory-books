function R = compare_coders(text)
%COMPARE_CODERS  Payload and header bits of the four coders, with
%   round-trip checks (header: 16 bits + 24 bits per symbol, as in Python).
[sym, ~, idx] = unique(text); counts = accumarray(idx(:), 1)'; p = counts / numel(text);
hdr = 16 + numel(sym) * 24;
sf = shannon_fano_code(sym, p);  b = prefix_encode(text, sym, sf);
assert(strcmp(prefix_decode(b, sym, sf), text)); R.sf = [numel(b) hdr];
hf = huffman_code(sym, p);       b = prefix_encode(text, sym, hf);
assert(strcmp(prefix_decode(b, sym, hf), text)); R.huffman = [numel(b) hdr];
b = arith_encode(text, sym, counts);
assert(strcmp(arith_decode(b, sym, counts, numel(text)), text)); R.arith = [numel(b) hdr];
c = lzw_encode(text, '');
assert(strcmp(lzw_decode(c, ''), text)); R.lzw = [lzw_bits(c, 256) 0];
R.H = entropy_bits(p); R.n = numel(text);
end
