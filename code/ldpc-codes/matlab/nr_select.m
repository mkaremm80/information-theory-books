function [bg, Z, ils] = nr_select(K, R, bg)
%NR_SELECT Base graph, lifting size and set index for K bits (data + CRC) at
% code rate R (TS 38.212, Sections 7.2.2 and 5.2.2).
if nargin < 3 || isempty(bg)
    if K <= 292 || (K <= 3824 && R <= 0.67) || R <= 0.25, bg = 2; else, bg = 1; end
end
if bg == 1
    kb = 22;
elseif K > 640, kb = 10; elseif K > 560, kb = 9; elseif K > 192, kb = 8; else, kb = 6;
end
S = nr_lifting_sets(); Z = inf; ils = -1;
for i = 1:8
    z = S{i}(find(kb * S{i} >= K, 1));
    if ~isempty(z) && z < Z, Z = z; ils = i - 1; end
end
end
