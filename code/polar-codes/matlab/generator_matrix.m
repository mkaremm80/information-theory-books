function G = generator_matrix(N)
%GENERATOR_MATRIX  G_N = F^(kron n) with F = [1 0; 1 1]  (for testing).
F = [1 0; 1 1];
G = 1;
while size(G,1) < N
    G = mod(kron(F, G), 2);
end
end
