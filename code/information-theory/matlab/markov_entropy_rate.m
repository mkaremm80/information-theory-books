function H = markov_entropy_rate(P)
%MARKOV_ENTROPY_RATE  H = sum_i pi_i * H(row i), in bits/symbol.
pi_ = markov_stationary(P);
Hr = zeros(1, size(P, 1));
for i = 1:size(P, 1)
    row = P(i, P(i, :) > 0);
    Hr(i) = -sum(row .* log2(row));
end
H = sum(pi_ .* Hr);
end
