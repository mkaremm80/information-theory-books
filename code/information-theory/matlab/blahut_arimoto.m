function [C, px] = blahut_arimoto(W, iters)
%BLAHUT_ARIMOTO  Capacity of a DMC W(x,y) = P(y|x) and the best input.
if nargin < 2, iters = 2000; end
nx = size(W, 1); px = ones(1, nx) / nx;
for it = 1:iters
    py = px * W;
    D = zeros(1, nx);
    for x = 1:nx
        k = W(x, :) > 0;
        D(x) = sum(W(x, k) .* log2(W(x, k) ./ py(k)));
    end
    new = px .* 2.^D; new = new / sum(new);
    if max(abs(new - px)) < 1e-12, px = new; break; end
    px = new;
end
C = channel_mi(px, W);
end
