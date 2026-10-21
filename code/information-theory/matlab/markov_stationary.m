function pi_ = markov_stationary(P)
%MARKOV_STATIONARY  Solve pi*P = pi, sum(pi) = 1.
S = size(P, 1);
A = [P' - eye(S); ones(1, S)];
pi_ = (A \ [zeros(S, 1); 1])';
end
