function s = sigma_from_ebn0(ebn0_db, R)
%SIGMA_FROM_EBN0 Noise standard deviation for BPSK (Es = 1).
s = sqrt(1 ./ (2 * R * 10 .^ (ebn0_db / 10)));
end
