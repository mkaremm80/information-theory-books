function codes = shannon_fano_code(sym, p)
%SHANNON_FANO_CODE  Shannon-Fano codes. sym: char array, p: probabilities.
%   Returns a cell array of codeword strings in the order of sym.
%   Sort by decreasing p (ties by symbol); split where |P_left - P_right|
%   is smallest (first such split); left part gets '0', right part '1'.
[~, order] = sortrows([-p(:), double(sym(:))]);
codes = repmat({''}, 1, numel(sym));
codes = split_group(order(:)', p, codes);
end

function codes = split_group(g, p, codes)
if numel(g) < 2, return; end
total = sum(p(g)); run = 0; best = inf; k = 1;
for j = 1:numel(g) - 1
    run = run + p(g(j));
    d = abs(2*run - total);
    if d < best - 1e-12, best = d; k = j; end
end
for j = 1:k,          codes{g(j)} = [codes{g(j)} '0']; end
for j = k+1:numel(g), codes{g(j)} = [codes{g(j)} '1']; end
codes = split_group(g(1:k), p, codes);
codes = split_group(g(k+1:end), p, codes);
end
