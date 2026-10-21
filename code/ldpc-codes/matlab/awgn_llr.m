function L = awgn_llr(c, ebn0_db, R)
%AWGN_LLR BPSK (0 -> +1, 1 -> -1) over AWGN; returns channel LLRs 2y/sigma^2.
sigma = sigma_from_ebn0(ebn0_db, R);
y = (1 - 2 * double(c)) + sigma * randn(size(c));
L = 2 * y / sigma ^ 2;
end
