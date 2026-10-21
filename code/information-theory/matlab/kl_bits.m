function D = kl_bits(p, q)
%KL_BITS  Kullback-Leibler divergence D(p||q) in bits (Inf if q=0 where p>0).
p = p(:); q = q(:);
k = p > 0;
if any(q(k) == 0), D = Inf; return; end
D = sum(p(k) .* log2(p(k) ./ q(k)));
end
