function z = bec_bhattacharyya(N, eps)
%BEC_BHATTACHARYYA  Exact Bhattacharyya parameters (erasure probabilities)
%   of the N bit-channels of a BEC(eps), natural index order (index i -> z(i+1)).
%       minus (worse)  channel:  Z -> 2Z - Z^2
%       plus  (better) channel:  Z -> Z^2
if N == 1
    z = eps;
    return
end
z = [bec_bhattacharyya(N/2, 2*eps - eps^2), bec_bhattacharyya(N/2, eps^2)];
end
