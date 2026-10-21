function c = bch1517_encode(m)
%BCH1517_ENCODE  Systematic BCH(15,7): c = [parity(8) message(7)],
%   coefficient order x^0..x^14, g(x) = 1 + x^4 + x^6 + x^7 + x^8.
g = [1 0 0 0 1 0 1 1 1];
r = [zeros(1, 8) m(:)'];                       % x^8 m(x)
for i = 15:-1:9                                % long division
    if r(i)
        r(i-8:i) = xor(r(i-8:i), g);
    end
end
c = [double(r(1:8)) m(:)'];
end
