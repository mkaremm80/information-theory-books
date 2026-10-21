function r = info_measures(P)
%INFO_MEASURES  Joint table P(x,y) (rows = X, columns = Y) -> struct with
%   HX, HY, HXY, HY_given_X, HX_given_Y, I (all in bits).
px = sum(P, 2);  py = sum(P, 1);
r.HX  = entropy_bits(px);
r.HY  = entropy_bits(py);
r.HXY = entropy_bits(P(:));
r.HY_given_X = r.HXY - r.HX;
r.HX_given_Y = r.HXY - r.HY;
r.I = r.HX + r.HY - r.HXY;
end
