function msg = prefix_decode(bits, sym, codes)
%PREFIX_DECODE  Decode a prefix-free code by reading bits until a match.
msg = ''; cur = '';
for b = bits
    cur = [cur b]; %#ok<AGROW>
    k = find(strcmp(codes, cur), 1);
    if ~isempty(k)
        msg = [msg sym(k)]; %#ok<AGROW>
        cur = '';
    end
end
if ~isempty(cur), error('bitstream ended inside a codeword'); end
end
