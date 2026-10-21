function c = encode_G(G, u)
%ENCODE_G Encoding with a generator matrix: c = u*G (mod 2).
c = mod(double(u) * G, 2);
end
