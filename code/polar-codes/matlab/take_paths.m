function A = take_paths(A, p)
%TAKE_PATHS  Re-order the paths (2nd dimension) of a B-by-L-by-n array:
%   A(b, l, :) <- A(b, p(b,l), :)   for every block b and path l.
[B, L, n] = size(A);
base = repmat((1:B).', 1, L) + (p - 1) * B;          % B-by-L linear index
idx = base + reshape((0:n-1) * B * L, 1, 1, n);      % implicit expansion
A = A(idx);
end
