function code = qc_lift(base, Z)
%QC_LIFT Expand a base matrix: -1 -> zero block, s >= 0 -> identity shifted
% right by s (row r has its one in column (r+s) mod Z, 0-based inside the block).
[I, J] = find(base >= 0);
r = (0:Z - 1).';
cn = zeros(Z * numel(I), 1); vn = cn;
for t = 1:numel(I)
    s = base(I(t), J(t));
    cn((t - 1) * Z + (1:Z)) = (I(t) - 1) * Z + r + 1;
    vn((t - 1) * Z + (1:Z)) = (J(t) - 1) * Z + mod(r + s, Z) + 1;
end
code = ldpc_from_edges(size(base, 1) * Z, size(base, 2) * Z, cn, vn, sprintf('QC(Z=%d)', Z));
end
