%FINAL_PROJECT_LAB  Information Theory and Compression Laboratory (MATLAB
%   version of project.py): statistics, four compressors, channel coding,
%   BSC errors, recovery.
text = repmat(['INFORMATION THEORY MEASURES SURPRISE. COMPRESSION REMOVES REDUNDANCY. ' ...
               'ERROR CORRECTION ADDS REDUNDANCY BACK, IN A CONTROLLED WAY. '], 1, 4);
[sym, ~, idx] = unique(text); counts = accumarray(idx(:), 1)'; p = counts / numel(text);
fprintf('STEP 1-3  %d symbols, %d distinct, H = %.4f bits/symbol\n', numel(text), numel(sym), entropy_bits(p));
R = compare_coders(text);
fprintf('STEP 4-8  SF %d+%d, Huffman %d+%d, arithmetic %d+%d, LZW %d (raw %d)\n', ...
        R.sf, R.huffman, R.arith, R.lzw(1), 8*numel(text));
hf = huffman_code(sym, p);
src = prefix_encode(text, sym, hf) - '0';
pBSC = 0.01; rng(2026);
schemes = {'none', 'hamming74', 'rep5'};
for s = 1:3
    switch schemes{s}
        case 'none',      tx = src;
        case 'hamming74', pad = [src zeros(1, mod(-numel(src), 4))];
                          tx = cell2mat(arrayfun(@(i) hamming74_encode(pad(i:i+3)), 1:4:numel(pad), 'UniformOutput', false));
        case 'rep5',      tx = repelem(src, 5);
    end
    rx = double(xor(tx, rand(size(tx)) < pBSC));
    switch schemes{s}
        case 'none',      out = rx;
        case 'hamming74', out = cell2mat(arrayfun(@(i) hamming74_decode(rx(i:i+6)), 1:7:numel(rx), 'UniformOutput', false));
        case 'rep5',      out = rep_decode(rx);
    end
    out = out(1:numel(src));
    ber = mean(out ~= src);
    try
        back = prefix_decode(char(out + '0'), sym, hf);
    catch
        back = '';                                 % broken tail after errors
    end
    fprintf('STEP 9-12 %-9s sent %5d bits, residual BER %.4f, text recovered exactly: %d\n', ...
            schemes{s}, numel(tx), ber, strcmp(back, text));
end
