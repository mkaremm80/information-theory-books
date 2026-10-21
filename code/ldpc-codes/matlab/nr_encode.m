function [c, x] = nr_encode(nr, u)
%NR_ENCODE Encode messages u (B x K): full codewords c and the E sent bits x.
s = zeros(size(u, 1), nr.Kp);
s(:, 1:nr.K) = u;                                      % filler bits stay 0
c = dual_encode(nr.enc, s);
x = c(:, nr.tx_pos);
end
