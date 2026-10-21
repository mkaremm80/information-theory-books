%LAB10_CHANNEL  Channel capacity (MATLAB version of channel.py).
h2 = @(p) -p.*log2(max(p, eps)) - (1-p).*log2(max(1-p, eps));
for p = [0 0.01 0.11 0.5]
    fprintf('BSC p = %.2f  C = %.4f bits/use\n', p, 1 - h2(p));
end
W = [0.9 0.1; 0.2 0.8];
[C, px] = blahut_arimoto(W);
fprintf('channel [0.9 0.1; 0.2 0.8]: C = %.4f, best P(X) = %s\n', C, mat2str(round(px*1e4)/1e4));
[C, px] = blahut_arimoto([1 0; 0.5 0.5]);
fprintf('Z-channel q = 0.5: C = %.4f, best P(X) = %s\n', C, mat2str(round(px*1e4)/1e4));
