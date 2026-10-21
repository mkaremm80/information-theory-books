function code = make_layers(code, Z)
%MAKE_LAYERS Groups of checks that share no variable. With Z: every block row of
% Z checks is one layer (QC codes); otherwise greedy colouring (same as Python).
first = [1; cumsum(code.dc(1:end - 1)) + 1];
if nargin > 1
    groups = arrayfun(@(r) (r:r + Z - 1).', 1:Z:code.m, 'UniformOutput', false);
else
    groups = {}; used = {};
    for c = 1:code.m
        vs = code.vn(first(c):first(c) + code.dc(c) - 1);
        placed = false;
        for g = 1:numel(groups)
            if ~any(used{g}(vs))
                groups{g}(end + 1, 1) = c; used{g}(vs) = true; placed = true; break;
            end
        end
        if ~placed
            groups{end + 1} = c; u = false(code.n, 1); u(vs) = true; used{end + 1} = u; %#ok<AGROW>
        end
    end
end
code.layers = cell(numel(groups), 1);
for g = 1:numel(groups)
    cs = groups{g}; e = [];
    for c = cs.', e = [e; (first(c):first(c) + code.dc(c) - 1).']; end %#ok<AGROW>
    L.e = e; L.m = numel(cs); L.cn = repelem((1:L.m).', code.dc(cs));
    L.Sc = sparse(1:numel(e), L.cn, 1, numel(e), L.m);
    L.dmax = max(code.dc(cs));
    t = (1:numel(e)).' - [0; cumsum(code.dc(cs(1:end - 1)))](L.cn);
    L.epos = L.cn + (t - 1) * L.m;
    code.layers{g} = L;
end
end
