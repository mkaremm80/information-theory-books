%LAB2_LAB3_CODES  Shannon-Fano and Huffman (MATLAB versions of Labs 2 and 3).
sym = 'ABCDEF';  p = [0.30 0.25 0.15 0.12 0.10 0.08];
H = entropy_bits(p);
sf = shannon_fano_code(sym, p);
[hf, merges] = huffman_code(sym, p);
for k = 1:numel(sym)
    fprintf('%s  p = %.2f   Shannon-Fano %-4s  Huffman %-4s\n', sym(k), p(k), sf{k}, hf{k});
end
Lsf = sum(p .* cellfun(@numel, sf));  Lhf = sum(p .* cellfun(@numel, hf));
fprintf('H = %.4f  L(SF) = %.4f  L(Huffman) = %.4f  efficiency %.2f%%\n', H, Lsf, Lhf, 100*H/Lhf);
for k = 1:size(merges, 1)
    fprintf('merge %d: %s(%.2f) + %s(%.2f)\n', k, merges{k,1}, merges{k,2}, merges{k,3}, merges{k,4});
end
msg = 'ABACADEF';
b = prefix_encode(msg, sym, hf);
fprintf('Huffman: %s -> %s -> %s\n', msg, b, prefix_decode(b, sym, hf));
