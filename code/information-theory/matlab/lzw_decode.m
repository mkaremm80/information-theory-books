function text = lzw_decode(codes, alphabet)
%LZW_DECODE  Rebuilds the dictionary while decoding; handles the cScSc case.
if nargin < 2 || isempty(alphabet), alphabet = char(0:255); end
d = cell(1, numel(alphabet) + numel(codes));
for i = 1:numel(alphabet), d{i} = alphabet(i); end
n = numel(alphabet);                         % number of entries
w = d{codes(1) + 1}; text = w;
for k = codes(2:end)
    if k < n
        entry = d{k + 1};
    elseif k == n
        entry = [w w(1)];                    % special case
    else
        error('bad code %d', k);
    end
    n = n + 1; d{n} = [w entry(1)];
    text = [text entry]; %#ok<AGROW>
    w = entry;
end
end
