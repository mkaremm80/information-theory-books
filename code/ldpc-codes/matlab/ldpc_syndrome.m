function s = ldpc_syndrome(code, C)
%LDPC_SYNDROME Syndromes H*c' for a batch of words C (B x n) -> (B x m).
s = mod(double(C(:, code.vn)) * code.Sc, 2);
end
