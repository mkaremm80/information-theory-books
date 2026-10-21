%LAB9_MARKOV  Markov chains (MATLAB version of markov.py).
P = [0.9 0.1; 0.5 0.5];
pi_ = markov_stationary(P);
fprintf('stationary distribution: %s\n', mat2str(round(pi_*1e4)/1e4));
fprintf('entropy rate: %.4f bits/symbol\n', markov_entropy_rate(P));
fprintf('entropy of pi (memoryless view): %.4f\n', entropy_bits(pi_));
rng(1); seq = markov_simulate(P, 20000, 1);
fprintf('estimated P from 20000 samples: %s\n', mat2str(round(markov_estimate(seq, 2)*1e3)/1e3));
