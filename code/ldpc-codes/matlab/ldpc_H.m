function H = ldpc_H(code)
%LDPC_H Dense parity-check matrix (small codes only).
H = full(sparse(code.cn, code.vn, 1, code.m, code.n));
end
