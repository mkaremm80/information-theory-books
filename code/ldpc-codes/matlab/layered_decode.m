function [chat, iters, app] = layered_decode(code, llr, max_iter, rule, alpha, beta)
%LAYERED_DECODE Row-layered decoder: the a-posteriori LLRs are updated after
% every layer, so new information is used at once. Call make_layers first.
if nargin < 4, rule = 'nms'; end
if nargin < 5, alpha = 0.75; end
if nargin < 6, beta = 0.5; end
B = size(llr, 1); app = llr; iters = max_iter * ones(B, 1);
act = (1:B).'; a = app; c2v = zeros(B, code.E);
for it = 1:max_iter
    for g = 1:numel(code.layers)
        L = code.layers{g}; e = L.e;
        v2c = a(:, code.vn(e)) - c2v(:, e);            % remove the old message
        nw = check_update(L, v2c, rule, alpha, beta);
        a(:, code.vn(e)) = v2c + nw;                   % update APP immediately
        c2v(:, e) = nw;
    end
    ok = ~any(ldpc_syndrome(code, a < 0), 2);
    app(act, :) = a;
    iters(act(ok)) = it;
    act = act(~ok); a = a(~ok, :); c2v = c2v(~ok, :);
    if isempty(act), break; end
end
chat = double(app < 0);
end
