function H = entropy_bits(p)
%ENTROPY_BITS  Shannon entropy H = -sum p log2 p (bits); zero probabilities skipped.
p = p(:);
if abs(sum(p) - 1) > 1e-9, error('probabilities must sum to 1'); end
p = p(p > 0);
H = -sum(p .* log2(p));
end
