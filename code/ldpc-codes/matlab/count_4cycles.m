function c = count_4cycles(code)
%COUNT_4CYCLES Number of 4-cycles = pairs of checks sharing two or more variables.
H = sparse(code.cn, code.vn, 1, code.m, code.n);
O = triu(full(H * H.'), 1);
c = sum(O(:) .* (O(:) - 1) / 2);
end
