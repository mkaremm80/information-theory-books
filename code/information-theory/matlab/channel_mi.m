function I = channel_mi(px, W)
%CHANNEL_MI  I(X;Y) for input px (row) and channel W(x,y) = P(y|x).
py = px * W;
I = 0;
for x = 1:numel(px)
    for y = 1:size(W, 2)
        if px(x) > 0 && W(x, y) > 0
            I = I + px(x) * W(x, y) * log2(W(x, y) / py(y));
        end
    end
end
end
