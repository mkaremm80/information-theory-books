function s = threshold_ga(lam, rho)
%THRESHOLD_GA Gaussian-approximation threshold (noise sigma) on BI-AWGN.
mt = logspace(-8, 3, 6000).';                          % table for the inverse
lp = log(max(phi_ga(mt), 1e-300));
[xl, ia] = unique(-lp);                               % increasing, no repeated values
lm = log(mt(ia));
inv = @(y) exp(interp1(xl, lm, -log(max(y, 1e-300)), 'linear', 'extrap'));
lo = 0.3; hi = 1.5;
for k = 1:30
    mid = (lo + hi) / 2;
    if de_ga(lam, rho, mid, inv), lo = mid; else, hi = mid; end
end
s = lo;
end

function ok = de_ga(lam, rho, sigma, inv)
m0 = 2 / sigma ^ 2; mu = 0; ok = false;
for it = 1:20000
    s = sum(lam(:, 2) .* phi_ga(m0 + (lam(:, 1) - 1) * mu));
    mu = 0;
    for r = 1:size(rho, 1)
        mu = mu + rho(r, 2) * inv(1 - (1 - s) ^ (rho(r, 1) - 1));
    end
    if mu > 60, ok = true; return; end
end
end
