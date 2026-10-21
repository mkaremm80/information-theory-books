%LAB1_ENTROPY  Entropy calculator (MATLAB version of entropy.py).
text = 'ABRACADABRA';
[sym, ~, idx] = unique(text);
p = accumarray(idx(:), 1)' / numel(text);
fprintf('symbols       : %s\n', sym);
fprintf('probabilities : %s\n', mat2str(round(p*1e4)/1e4));
fprintf('H(X)          : %.4f bits/symbol\n', entropy_bits(p));
fprintf('max possible  : %.4f bits/symbol\n', log2(numel(sym)));
fprintf('H2(0.1)       : %.4f bits\n', entropy_bits([0.1 0.9]));
r = info_measures([0.25 0.25; 0 0.5]);
fprintf('H(X,Y) = %.4f, H(Y|X) = %.4f, I(X;Y) = %.4f\n', r.HXY, r.HY_given_X, r.I);
fprintf('min-entropy of [0.6 0.2 0.1 0.1]: %.4f bits\n', -log2(0.6));
