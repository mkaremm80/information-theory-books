function [low, high, tag] = arith_interval(msg, sym, p)
%ARITH_INTERVAL  Teaching version: interval of a short message and the
%   shortest binary tag inside it (double precision: use for < 40 symbols).
[sym, o] = sort(sym); p = p(o);
Fh = cumsum(p); Fl = Fh - p;
low = 0; high = 1;
for c = msg
    k = find(sym == c, 1);
    w = high - low;
    high = low + w * Fh(k);
    low  = low + w * Fl(k);
    fprintf('  after %s: [%.6f, %.6f)  width %.6f\n', c, low, high, high - low);
end
k = 1;
while true
    m = ceil(low * 2^k - 1e-12);
    if (m + 1) / 2^k <= high + 1e-12
        tag = dec2bin(m, k); return
    end
    k = k + 1;
end
end
