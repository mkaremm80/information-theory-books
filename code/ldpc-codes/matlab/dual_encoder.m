function enc = dual_encoder(code, k, g)
%DUAL_ENCODER Prepare the encoder for H = [A B 0; C1 C2 I]:
% k systematic bits, g core parity bits, m-g extension parity bits.
enc.code = code; enc.k = k; enc.g = g;
sel = code.cn <= g & code.vn <= k + g;
Hc = full(sparse(code.cn(sel), code.vn(sel), 1, g, k + g));
enc.Binv = gf2_inv(Hc(:, k + 1:k + g));                % small g x g inverse
diagE = code.cn > g & code.vn == code.cn + k;          % the identity part
enc.ext = code.cn > g & ~diagE;                        % edges of [C1 C2]
enc.core = code.cn <= g & code.vn <= k;                % edges of A
end
