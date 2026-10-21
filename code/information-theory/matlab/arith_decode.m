function msg = arith_decode(bits, sym, counts, n)
%ARITH_DECODE  Decoder for ARITH_ENCODE; n = number of symbols to decode.
[sym, o] = sort(sym); counts = counts(o);
cumh = cumsum(counts); cuml = cumh - counts; total = cumh(end);
FULL = 2^32 - 1; HALF = 2^31; QUART = 2^30;
stream = bits - '0'; pos = 0;
low = 0; high = FULL; code = 0;
for i = 1:32, [b, pos] = nextbit(stream, pos); code = 2*code + b; end
msg = blanks(n);
for i = 1:n
    rng = high - low + 1;
    target = floor(((code - low + 1) * total - 1) / rng);
    k = find(target >= cuml & target < cumh, 1);
    msg(i) = sym(k);
    high = low + floor(rng * cumh(k) / total) - 1;
    low  = low + floor(rng * cuml(k) / total);
    while true
        if high < HALF
            % nothing to subtract
        elseif low >= HALF
            low = low - HALF; high = high - HALF; code = code - HALF;
        elseif low >= QUART && high < 3*QUART
            low = low - QUART; high = high - QUART; code = code - QUART;
        else
            break
        end
        low = 2*low; high = 2*high + 1;
        [b, pos] = nextbit(stream, pos); code = 2*code + b;
    end
end
end

function [b, pos] = nextbit(stream, pos)
pos = pos + 1;
if pos <= numel(stream), b = stream(pos); else, b = 0; end
end
