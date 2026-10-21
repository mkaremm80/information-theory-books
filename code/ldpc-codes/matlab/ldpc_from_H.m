function code = ldpc_from_H(H, name)
%LDPC_FROM_H Build the Tanner-graph object from a dense parity-check matrix.
if nargin < 2, name = 'LDPC'; end
[cn, vn] = find(mod(H, 2));
code = ldpc_from_edges(size(H, 1), size(H, 2), cn, vn, name);
end
