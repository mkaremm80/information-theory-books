function code = ldpc_from_edges(m, n, cn, vn, name)
%LDPC_FROM_EDGES Tanner-graph object. cn, vn: check and variable of each edge
% (1-based). Edges are sorted by check, then by variable (same order as Python).
if nargin < 5, name = 'LDPC'; end
cv = sortrows([cn(:), vn(:)]);
code.m = m; code.n = n; code.name = name;
code.cn = cv(:, 1); code.vn = cv(:, 2); code.E = numel(code.cn);
code.dc = accumarray(code.cn, 1, [m 1]);           % check degrees
code.dv = accumarray(code.vn, 1, [n 1]);           % variable degrees
code.Sc = sparse(1:code.E, code.cn, 1, code.E, m); % edge -> check incidence
code.Sv = sparse(1:code.E, code.vn, 1, code.E, n); % edge -> variable incidence
code = add_padding(code);
code.layers = {};
end

function s = add_padding(s)
% position of every edge inside an (m x max degree) table, for min-sum
first = [1; cumsum(s.dc(1:end - 1)) + 1];
t = (1:s.E).' - first(s.cn) + 1;                   % slot of the edge in its check
s.dmax = max(s.dc);
s.epos = s.cn + (t - 1) * s.m;                     % linear index in m x dmax
end
