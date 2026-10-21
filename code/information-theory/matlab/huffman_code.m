function [codes, merges] = huffman_code(sym, p)
%HUFFMAN_CODE  Huffman codes with the book's conventions:
%   nodes are kept with an 'age'; the two smallest (ties: older first) are
%   merged, the first removed gets bit '0', the second bit '1'.
n = numel(sym);
codes = repmat({''}, 1, n);
if n == 1, codes{1} = '0'; merges = {}; return; end
[~, alpha] = sort(double(sym));                 % leaves enter in symbol order
prob = p(alpha); members = num2cell(alpha); age = 1:n; nextage = n + 1;
merges = {};
while numel(prob) > 1
    [~, o] = sortrows([prob(:), age(:)]);         % smallest, then oldest
    a = o(1); b = o(2);
    for s = members{a}, codes{s} = ['0' codes{s}]; end
    for s = members{b}, codes{s} = ['1' codes{s}]; end
    merges(end+1, :) = {sym(members{a}), prob(a), sym(members{b}), prob(b)}; %#ok<AGROW>
    newp = prob(a) + prob(b); newm = [members{a} members{b}];
    keep = setdiff(1:numel(prob), [a b]);
    prob = [prob(keep) newp]; members = [members(keep) {newm}];
    age = [age(keep) nextage]; nextage = nextage + 1;
end
end
