function y = phi_ga(m)
%PHI_GA phi(m) = 1 - E[tanh(u/2)], u ~ N(m, 2m) (Chung et al. 2001), by quadrature.
persistent x w
if isempty(x), [x, w] = gh_nodes(80); end
m = m(:);
u = m + sqrt(2 * m) .* x.';
y = 1 - (tanh(u / 2) * w) / sqrt(2 * pi);
y(m < 1e-10) = 1;
end
