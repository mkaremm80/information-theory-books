function [u_hat, u_list, PM] = polar_list_decode(llr, info, L, crc)
%POLAR_LIST_DECODE  LLR-based SC / SCL / CRC-aided SCL polar decoder.
%   llr  : B-by-N channel LLRs (one block per row), L = ln P(0)/P(1).
%   info : 1-by-N logical information set (frozen bits are 0).
%   L    : list size. L = 1 gives the plain SC decoder.
%   crc  : (optional) struct from crc_make(); enables CA-SCL selection.
%   u_hat  : B-by-N decided u vectors.
%   u_list : B-by-L-by-N final list, sorted by path metric (best first).
%   PM     : B-by-L path metrics (min-sum form, smaller = better).
%
%   Path metric update (Balatsoukas-Stimming et al., 2015):
%       PM <- PM + |LLR|  if the decision disagrees with sign(LLR).
if nargin < 3, L = 1; end
if nargin < 4, crc = []; end
[B, N] = size(llr);
st.info = info;
st.PM = inf(B, L);
st.PM(:, 1) = 0;                                  % one active path at start
lam = repmat(reshape(llr, B, 1, N), 1, L, 1);     % B-by-L-by-N
[~, u_list, ~, st] = node(lam, 1, st);
PM = st.PM;
if isempty(crc)
    u_hat = reshape(u_list(:, 1, :), B, N);       % best path metric
else
    payload = u_list(:, :, info);                 % B-by-L-by-K
    ok = reshape(crc_check(crc, reshape(payload, B*L, [])), B, L) & isfinite(PM);
    pick = ones(B, 1);
    for b = 1:B
        f = find(ok(b, :), 1);
        if ~isempty(f), pick(b) = f; end
    end
    idx = sub2ind([B, L], (1:B).', pick);
    U = reshape(u_list, B*L, N);
    u_hat = U(idx, :);
end
end

% ---------------------------------------------------------------- tree node
function [c, u, parent, st] = node(lam, start, st)
% lam: B-by-L-by-n LLRs of this node; start: index of its first leaf (1-based)
n = size(lam, 3);
if n == 1
    [c, u, parent, st] = leaf(lam, start, st);
    return
end
h = n/2;
a = lam(:, :, 1:h);  b = lam(:, :, h+1:n);
% left child: f (min-sum)
[c1, u1, p1, st] = node(sign(a).*sign(b).*min(abs(a), abs(b)), start, st);
% right child: g, with the surviving parents' LLRs
a1 = take_paths(a, p1);  b1 = take_paths(b, p1);
[c2, u2, p2, st] = node(b1 + (1 - 2*c1).*a1, start + h, st);
% re-order left results by the right child's survivors
c1 = take_paths(c1, p2);
u1 = take_paths(u1, p2);
parent = p1(sub2ind(size(p1), repmat((1:size(p1,1)).', 1, size(p1,2)), p2));
c = cat(3, double(xor(c1, c2)), c2);      % partial sums [c1 xor c2 | c2]
u = cat(3, u1, u2);
end

% ---------------------------------------------------------------- leaf
function [c, u, parent, st] = leaf(lam, i, st)
llr = lam(:, :, 1);                        % B-by-L
[B, L] = size(llr);
pen0 = (llr < 0) .* (-llr);                % cost of deciding 0
pen1 = (llr > 0) .* llr;                   % cost of deciding 1
if ~st.info(i)                             % frozen bit: always 0
    st.PM = st.PM + pen0;
    c = zeros(B, L);
    parent = repmat(1:L, B, 1);
else                                       % information bit: split and prune
    cand = [st.PM + pen0, st.PM + pen1];   % B-by-2L
    [srt, ord] = sort(cand, 2);            % stable sort
    ord = ord(:, 1:L);
    st.PM = srt(:, 1:L);
    parent = mod(ord - 1, L) + 1;
    c = double(ord > L);
end
u = c;
end
