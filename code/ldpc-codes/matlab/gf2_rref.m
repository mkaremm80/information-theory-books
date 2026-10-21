function [R, piv] = gf2_rref(A)
%GF2_RREF Reduced row echelon form over GF(2). piv = pivot columns (1-based).
R = mod(double(A), 2);
[m, n] = size(R);
piv = []; row = 1;
for col = 1:n
    if row > m, break; end
    hits = find(R(row:end, col), 1);
    if isempty(hits), continue; end
    p = row + hits - 1;
    R([row p], :) = R([p row], :);                 % swap the pivot row up
    others = find(R(:, col)); others(others == row) = [];
    R(others, :) = mod(R(others, :) + R(row, :), 2);   % clear the column
    piv(end + 1) = col; %#ok<AGROW>
    row = row + 1;
end
R = R(1:row - 1, :);
end
