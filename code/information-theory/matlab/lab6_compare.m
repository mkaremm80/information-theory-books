%LAB6_COMPARE  Fair comparison on the four data sets of Chapter 10
%   (the text files were written by compare.py so both languages use the same data).
names = {'skewed', 'english', 'repetitive', 'uniform'};
for k = 1:numel(names)
    fid = fopen(['data_' names{k} '.txt']); text = fgetl(fid); fclose(fid);
    R = compare_coders(text);
    fprintf('%-10s n = %4d  H = %.3f  nH = %4.0f |', names{k}, R.n, R.H, R.n * R.H);
    fprintf(' SF %d+%d | Huffman %d+%d | Arith %d+%d | LZW %d\n', R.sf, R.huffman, R.arith, R.lzw(1));
end
