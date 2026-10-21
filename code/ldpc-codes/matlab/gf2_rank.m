function r = gf2_rank(A)
%GF2_RANK Rank of a binary matrix over GF(2).
[~, piv] = gf2_rref(A);
r = numel(piv);
end
