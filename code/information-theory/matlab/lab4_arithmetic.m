%LAB4_ARITHMETIC  Arithmetic coding (MATLAB version of arithmetic.py).
disp('Exact example: p(A)=1/2, p(B)=1/4, p(C)=1/4, message ABAC');
[lo, hi, tag] = arith_interval('ABAC', 'ABC', [0.5 0.25 0.25]);
fprintf('  tag bits: %s\n', tag);
text = 'TOBEORNOTTOBEORTOBEORNOT';
[sym, ~, idx] = unique(text); counts = accumarray(idx(:), 1)';
b = arith_encode(text, sym, counts);
fprintf('integer coder: %d bits; round trip OK: %d\n', numel(b), ...
        strcmp(arith_decode(b, sym, counts, numel(text)), text));
