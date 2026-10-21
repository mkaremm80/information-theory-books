function bits = arith_encode(msg, sym, counts)
%ARITH_ENCODE  Integer arithmetic coder (32-bit, E1/E2/E3 scaling).
%   msg: char row; sym: char row of symbols (any order); counts: integer
%   frequencies (total <= 2^16 so that all products stay exact in doubles).
%   Symbols are placed in the cumulative table in sorted order.
[sym, o] = sort(sym); counts = counts(o);
cumh = cumsum(counts); cuml = cumh - counts; total = cumh(end);
FULL = 2^32 - 1; HALF = 2^31; QUART = 2^30;
low = 0; high = FULL; pending = 0; bits = blanks(0);
for c = msg
    k = find(sym == c, 1);
    rng = high - low + 1;
    high = low + floor(rng * cumh(k) / total) - 1;
    low  = low + floor(rng * cuml(k) / total);
    while true
        if high < HALF
            [bits, pending] = emit(bits, '0', pending);
        elseif low >= HALF
            [bits, pending] = emit(bits, '1', pending);
            low = low - HALF; high = high - HALF;
        elseif low >= QUART && high < 3*QUART
            pending = pending + 1; low = low - QUART; high = high - QUART;
        else
            break
        end
        low = 2*low; high = 2*high + 1;
    end
end
pending = pending + 1;
if low < QUART, [bits, ~] = emit(bits, '0', pending);
else,           [bits, ~] = emit(bits, '1', pending); end
end

function [bits, pending] = emit(bits, b, pending)
bits = [bits b repmat(char('0' + ('1' - b)), 1, pending)];
pending = 0;
end
