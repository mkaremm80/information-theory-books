function H = gallager_H(n, dv, dc)
%GALLAGER_H Gallager's (dv,dc)-regular construction: dv stacked bands, each a
% random column permutation of the first band (uses MATLAB's random stream).
mb = n / dc;
band = kron(eye(mb), ones(1, dc));
H = band;
for b = 2:dv
    H = [H; band(:, randperm(n))]; %#ok<AGROW>
end
end
