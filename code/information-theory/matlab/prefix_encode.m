function bits = prefix_encode(msg, sym, codes)
%PREFIX_ENCODE  Concatenate codewords of the symbols in msg.
bits = '';
for c = msg
    bits = [bits codes{sym == c}]; %#ok<AGROW>
end
end
