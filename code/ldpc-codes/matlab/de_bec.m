function hist = de_bec(lam, rho, eps, iters, tol)
%DE_BEC Density evolution on the BEC: x(l+1) = eps*lambda(1 - rho(1 - x(l))).
% lam, rho: [degree, edge fraction] rows (edge perspective).
if nargin < 4, iters = 2000; end
if nargin < 5, tol = 1e-12; end
x = eps; hist = eps;
for l = 1:iters
    x = eps * dpoly(lam, 1 - dpoly(rho, 1 - x));
    hist(end + 1) = x; %#ok<AGROW>
    if x < tol, break; end
end
end

function y = dpoly(d, x)
y = sum(d(:, 2) .* x .^ (d(:, 1) - 1));
end
