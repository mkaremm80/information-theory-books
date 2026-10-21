function [chat, iters, app] = bp_decode(code, llr, max_iter, rule, alpha, beta, early_stop)
%BP_DECODE Flooding belief propagation on a batch of frames (rows of llr).
% Returns hard decisions, iterations used per frame and a-posteriori LLRs.
if nargin < 4, rule = 'spa'; end
if nargin < 5, alpha = 0.75; end
if nargin < 6, beta = 0.5; end
if nargin < 7, early_stop = true; end
B = size(llr, 1);
app = llr; iters = max_iter * ones(B, 1);
act = (1:B).';                                         % frames still decoding
c2v = zeros(B, code.E);
for it = 1:max_iter
    v2c = app(act, code.vn) - c2v;                     % extrinsic: remove own message
    c2v = check_update(code, v2c, rule, alpha, beta);
    app(act, :) = llr(act, :) + c2v * code.Sv;         % a-posteriori LLRs
    if early_stop
        ok = ~any(ldpc_syndrome(code, app(act, :) < 0), 2);
        iters(act(ok)) = it;
        act = act(~ok); c2v = c2v(~ok, :);
        if isempty(act), break; end
    end
end
chat = double(app < 0);
end
