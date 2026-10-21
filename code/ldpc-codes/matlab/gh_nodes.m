function [x, w] = gh_nodes(N)
%GH_NODES Gauss-Hermite nodes and weights for the weight exp(-x^2/2)
% (probabilists' Hermite polynomials), by the Golub-Welsch method.
J = diag(sqrt(1:N - 1), 1) + diag(sqrt(1:N - 1), -1);
[V, D] = eig(J);
[x, i] = sort(diag(D));
w = sqrt(2 * pi) * V(1, i).' .^ 2;
end
