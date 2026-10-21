function [x, iters] = bit_flip(code, r, max_iter)
%BIT_FLIP Parallel bit flipping: flip every bit that sits in the largest number
% of failed checks. r: hard decisions (B x n).
x = double(r); iters = zeros(size(x, 1), 1);
for it = 1:max_iter
    s = ldpc_syndrome(code, x);
    active = any(s, 2);
    if ~any(active), break; end
    iters(active) = it;
    fails = s(:, code.cn) * code.Sv;                   % failed checks per bit
    worst = max(fails, [], 2);
    flip = (fails == worst) & (worst > 0) & active;
    x = mod(x + flip, 2);
end
end
