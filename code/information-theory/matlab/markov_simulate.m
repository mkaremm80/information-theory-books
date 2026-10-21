function seq = markov_simulate(P, n, start)
%MARKOV_SIMULATE  n states of a chain (states 1..S); seed rng before calling.
seq = zeros(1, n); seq(1) = start; C = cumsum(P, 2);
for k = 2:n
    seq(k) = find(rand < C(seq(k-1), :), 1);
end
end
