function t = threshold_bec(lam, rho, iters)
%THRESHOLD_BEC Largest erasure probability for which DE goes to zero (bisection).
if nargin < 3, iters = 5000; end
lo = 0; hi = 1;
for k = 1:40
    mid = (lo + hi) / 2;
    h = de_bec(lam, rho, mid, iters);
    if h(end) < 1e-10, lo = mid; else, hi = mid; end
end
t = lo;
end
