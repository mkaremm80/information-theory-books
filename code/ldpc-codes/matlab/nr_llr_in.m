function L = nr_llr_in(nr, llr_tx)
%NR_LLR_IN Decoder input: LLRs of the sent bits at their positions (repeated bits
% add up), 0 for punctured / unsent bits, a large value for filler bits.
L = zeros(size(llr_tx, 1), nr.n);
N0 = numel(nr.tx_pos0);
for p = 1:N0:nr.E
    q = min(p + N0 - 1, nr.E);
    L(:, nr.tx_pos(p:q)) = L(:, nr.tx_pos(p:q)) + llr_tx(:, p:q);
end
L(:, nr.filler) = 1e3;
end
