function Ai = gf2_inv(A)
%GF2_INV Inverse of a square invertible matrix over GF(2).
n = size(A, 1);
[R, piv] = gf2_rref([mod(A, 2), eye(n)]);
if numel(piv) < n || any(piv(1:n) ~= 1:n)
    error('matrix is singular over GF(2)');
end
Ai = R(:, n + 1:end);
end
