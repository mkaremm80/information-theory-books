function g = ldpc_girth(code)
%LDPC_GIRTH Shortest cycle of the Tanner graph (breadth-first search from
% every variable node). Nodes 1..n are variables, n+1..n+m are checks.
A = [sparse(code.n, code.n), sparse(code.vn, code.cn, 1, code.n, code.m); ...
     sparse(code.cn, code.vn, 1, code.m, code.n), sparse(code.m, code.m)];
nb = cell(code.n + code.m, 1);
for v = 1:code.n + code.m, nb{v} = find(A(:, v)).'; end
g = inf;
for root = 1:code.n
    dist = -ones(code.n + code.m, 1); par = zeros(code.n + code.m, 1);
    dist(root) = 0; queue = root; head = 1;
    while head <= numel(queue)
        u = queue(head); head = head + 1;
        if 2 * dist(u) >= g, break; end
        for w = nb{u}
            if w == par(u), continue; end
            if dist(w) >= 0
                g = min(g, dist(u) + dist(w) + 1);
            else
                dist(w) = dist(u) + 1; par(w) = u; queue(end + 1) = w; %#ok<AGROW>
            end
        end
    end
end
end
