function c = dual_encode(enc, s)
%DUAL_ENCODE Three-step encoding: lam = A s; p1 = B^-1 lam; p2 = C1 s + C2 p1.
code = enc.code; k = enc.k; g = enc.g;
c = zeros(size(s, 1), code.n);
c(:, 1:k) = s;
lam = rowsum(code, enc.core, c);                       % step 1
c(:, k + 1:k + g) = mod(lam(:, 1:g) * enc.Binv.', 2);  % step 2
p2 = rowsum(code, enc.ext, c);                         % step 3
c(:, k + g + 1:end) = p2(:, g + 1:code.m);
end

function out = rowsum(code, mask, c)
e = find(mask);
S = sparse(1:numel(e), code.cn(e), 1, numel(e), code.m);
out = mod(c(:, code.vn(e)) * S, 2);
end
