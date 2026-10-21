function [y, rounds] = peel(code, y)
%PEEL Peeling decoder for the BEC. y holds 0, 1 or -1 (erasure). A check with
% exactly one erased neighbour recovers it (all checks in parallel each round).
rounds = zeros(size(y, 1), 1);
for it = 1:code.n
    ev = y(:, code.vn);
    n_er = double(ev < 0) * code.Sc;                   % erasures per check
    par = mod(max(ev, 0) * code.Sc, 2);                % parity of known bits
    solve = n_er == 1;
    if ~any(solve(:)), break; end
    es = solve(:, code.cn) & (ev < 0);                 % the single erased edge
    [b, e] = find(es); b = b(:); e = e(:);
    y(sub2ind(size(y), b, code.vn(e))) = par(sub2ind(size(par), b, code.cn(e)));
    rounds(unique(b)) = it;
end
end
