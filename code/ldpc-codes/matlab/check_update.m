function out = check_update(seg, v2c, rule, alpha, beta)
%CHECK_UPDATE Check-node messages from variable-to-check messages v2c (B x E).
% seg describes the checks: fields cn (check of each edge), m, Sc, dmax, epos.
% rule: 'spa' (sum-product), 'ms', 'nms' (x alpha) or 'oms' (- beta).
if nargin < 4, alpha = 0.75; end
if nargin < 5, beta = 0.5; end
B = size(v2c, 1);
neg = double(v2c < 0);
sign_all = mod(neg * seg.Sc, 2);                      % parity of negatives per check
sgn = 1 - 2 * mod(sign_all(:, seg.cn) + neg, 2);      % exclude own sign
mag = abs(v2c);
if strcmp(rule, 'spa')
    p = phi_fn(mag);
    S = p * seg.Sc;                                    % sum of phi per check
    out = phi_fn(S(:, seg.cn) - p);
else
    m1 = groupmin(seg, mag, B);                        % smallest magnitude
    is1 = mag == m1(:, seg.cn);
    n1 = double(is1) * seg.Sc;                         % how many equal it
    mag2 = mag; mag2(is1) = inf;
    m2 = groupmin(seg, mag2, B);
    m2(n1 >= 2) = m1(n1 >= 2);                         % second smallest
    out = m1(:, seg.cn);
    M2 = m2(:, seg.cn);
    out(is1) = M2(is1);
    if strcmp(rule, 'nms')
        out = alpha * out;
    elseif strcmp(rule, 'oms')
        out = max(out - beta, 0);
    end
end
out = sgn .* out;
end

function mn = groupmin(seg, x, B)
T = inf(B, seg.m * seg.dmax);
T(:, seg.epos) = x;
mn = min(reshape(T, B, seg.m, seg.dmax), [], 3);
if B == 1, mn = mn(:).'; end
end
