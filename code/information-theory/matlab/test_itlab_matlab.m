%TEST_ITLAB_MATLAB  Tests of the MATLAB labs, including comparison with
%   reference results of the Python library (pyref_it.mat).
R = load('pyref_it.mat');
syms = 'ABCDEFGHIJ';
cellrow = @(x) strsplit(strtrim(x), ' ');
okSF = 0; okHF = 0;
for t = 1:50
    k = R.K(t); p = R.P(t, 1:k);
    sf = shannon_fano_code(syms(1:k), p);  hf = huffman_code(syms(1:k), p);
    okSF = okSF + isequal(sf, cellrow(R.SF{t}));
    okHF = okHF + isequal(hf, cellrow(R.HF{t}));
end
fprintf('Shannon-Fano codes identical to Python: %d / 50\n', okSF);
fprintf('Huffman codes identical to Python     : %d / 50\n', okHF);
okA = 0; okL = 0; okRT = 0;
for t = 1:30
    txt = R.T{t};
    [sym, ~, idx] = unique(txt); counts = accumarray(idx(:), 1)';
    b = arith_encode(txt, sym, counts);
    okA = okA + strcmp(b, R.AB{t});
    okRT = okRT + strcmp(arith_decode(b, sym, counts, numel(txt)), txt);
    c = lzw_encode(txt, 'ABCDE');
    okL = okL + isequal(c, str2num(R.LZ{t})) * strcmp(lzw_decode(c, 'ABCDE'), txt); %#ok<ST2NM>
end
fprintf('arithmetic bit streams identical      : %d / 30 (round trips %d / 30)\n', okA, okRT);
fprintf('LZW codes identical and round trip    : %d / 30\n', okL);
okB = 0; okE = 0;
for t = 1:20
    c = bch1517_encode(R.M(t, :));
    okB = okB + isequal(c, R.C(t, :));
    for i = 0:14
        for j = i:14
            r = c; r(i+1) = 1 - r(i+1);
            if j > i, r(j+1) = 1 - r(j+1); end
            [~, ch] = bch1517_decode(r);
            okE = okE + isequal(ch, c);
        end
    end
end
fprintf('BCH codewords identical               : %d / 20\n', okB);
fprintf('BCH 1- and 2-error patterns corrected : %d / %d\n', okE, 20 * 120);
okH = 0;
for v = 0:15
    d = bitget(v, 4:-1:1); c = hamming74_encode(d);
    for e = 1:7
        r = c; r(e) = 1 - r(e);
        okH = okH + isequal(hamming74_decode(r), d);
    end
end
fprintf('Hamming single errors corrected       : %d / 112\n', okH);
okR = 0;
for v = 0:31
    e = bitget(v, 1:5);
    okR = okR + (isequal(rep_decode(e), 0) == (sum(e) <= 2));
end
fprintf('(5,1) patterns behaving as predicted  : %d / 32\n', okR);
