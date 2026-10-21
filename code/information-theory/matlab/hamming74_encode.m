function c = hamming74_encode(d)
%HAMMING74_ENCODE  d = [d1 d2 d3 d4] -> [p1 p2 d1 p4 d2 d3 d4] (positions 1..7).
p1 = xor(xor(d(1), d(2)), d(4));
p2 = xor(xor(d(1), d(3)), d(4));
p4 = xor(xor(d(2), d(3)), d(4));
c = double([p1 p2 d(1) p4 d(2) d(3) d(4)]);
end
