function P = markov_estimate(seq, S, alpha)
%MARKOV_ESTIMATE  ML transition matrix from a state sequence (states 1..S),
%   with optional add-alpha smoothing.
if nargin < 3, alpha = 0; end
N = alpha * ones(S);
for k = 2:numel(seq)
    N(seq(k-1), seq(k)) = N(seq(k-1), seq(k)) + 1;
end
P = N ./ sum(N, 2);
end
