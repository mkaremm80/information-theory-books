function [G, info] = generator_from_H(H)
%GENERATOR_FROM_H Generator matrix G (k x n) with H*G' = 0 (mod 2).
% info = columns (1-based) where the codeword equals the message.
[R, piv] = gf2_rref(H);
n = size(H, 2);
info = setdiff(1:n, piv);
G = zeros(numel(info), n);
G(:, info) = eye(numel(info));
G(:, piv) = R(:, info).';            % c(piv_i) = sum_j R(i, info_j) c(info_j)
end
