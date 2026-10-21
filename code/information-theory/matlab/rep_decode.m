function [bits, votes] = rep_decode(r, n)
%REP_DECODE  Majority vote on blocks of n bits (n odd, default 5).
if nargin < 2, n = 5; end
blocks = reshape(r, n, []);
votes = sum(blocks, 1);
bits = double(votes > floor(n/2));
end
