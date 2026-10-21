function [info, m] = ga_construct(N, K, design_ebn0_db, rate)
%GA_CONSTRUCT  Gaussian-approximation construction for BPSK over AWGN.
%   info : 1-by-N logical, true at the K most reliable positions.
%   m    : mean LLR of every bit-channel (larger = more reliable).
%   rate : rate used to convert Eb/N0 to Es/N0 (default K/N).
%       plus  child:  m+ = 2m
%       minus child:  m- = phi^-1( 1 - (1 - phi(m))^2 )
if nargin < 3, design_ebn0_db = 2.0; end
if nargin < 4, rate = K/N; end
esn0 = 10^(design_ebn0_db/10) * rate;
m = ga_rec(N, 4*esn0);                   % channel LLR mean = 2/sigma^2
[~, order] = sort(m);
info = false(1, N);
info(order(end-K+1:end)) = true;
end

function m = ga_rec(N, m0)
% Level by level: each mean m becomes the pair (m-, m+). Interleaving the
% pairs gives natural index order (first split = most significant bit).
m = m0;
while numel(m) < N
    m = reshape([ga_minus(m); 2*m], 1, []);
end
end

function mm = ga_minus(m)
% Solve log(phi(mm)) = log(2p - p^2), p = phi(m), by bisection in the log
% domain (vectorised over all means of one level).
lp = log_phi(m);
p = exp(lp);
target = lp + log(2 - p);
lo = zeros(size(m)); hi = m;
for it = 1:80
    mid = (lo + hi)/2;
    up = log_phi(mid) > target;
    lo(up) = mid(up);
    hi(~up) = mid(~up);
end
mm = (lo + hi)/2;
end

function lp = log_phi(x)
% log of phi(x) = 1 - E[tanh(U/2)], U ~ N(x, 2x)   (vectorised).
%   x < 10 : exact value from a table built once by numerical integration
%   x >= 10: asymptotic form  sqrt(pi/x) exp(-x/4) (1 - 10/(7x))
persistent lxg lt
if isempty(lxg)
    xg = logspace(-8, 1, 4000);
    t = linspace(-12, 12, 4001);          % standard normal grid
    w = exp(-t.^2/2) / sqrt(2*pi);
    lt = zeros(size(xg));
    for k = 1:numel(xg)
        u = xg(k) + sqrt(2*xg(k)) * t;
        lt(k) = log(max(1 - trapz(t, w .* tanh(u/2)), 1e-300));
    end
    lxg = log(xg);
end
lp = zeros(size(x));
s = x <= 1e-8;            lp(s) = log1p(-x(s)/2);
s = x > 1e-8 & x < 10;    lp(s) = interp1(lxg, lt, log(x(s)));
s = x >= 10;              lp(s) = 0.5*log(pi./x(s)) - x(s)/4 + log(1 - 10./(7*x(s)));
end
