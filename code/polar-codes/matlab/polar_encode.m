function x = polar_encode(u)
%POLAR_ENCODE  Butterfly polar encoder, x = u*G_N (mod 2), G_N = F^(kron n).
%   u : 1-by-N vector or B-by-N matrix (one block per row) of 0/1 values.
%   x : same size as u.
%   Convention (as in 3GPP NR): F = [1 0; 1 1], natural order, no bit reversal.
%   Each stage XORs the lower half of every block into the upper half.
x = logical(u);
[B, N] = size(x);
h = N/2;
while h >= 1
    x = reshape(x.', 2*h, []);            % each column = one block of 2h bits
    x(1:h, :) = xor(x(1:h, :), x(h+1:2*h, :));
    x = reshape(x, N, B).';
    h = h/2;
end
x = double(x);
end
