function codes = lzw_encode(text, alphabet)
%LZW_ENCODE  LZW codes. alphabet = '' uses the 256 byte values (code = ASCII).
if nargin < 2 || isempty(alphabet), alphabet = char(0:255); end
d = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(alphabet), d(alphabet(i)) = i - 1; end
w = ''; codes = [];
for c = text
    if isKey(d, [w c])
        w = [w c];
    else
        codes(end+1) = d(w); %#ok<AGROW>
        d([w c]) = d.Count;      % next free code
        w = c;
    end
end
if ~isempty(w), codes(end+1) = d(w); end
end
