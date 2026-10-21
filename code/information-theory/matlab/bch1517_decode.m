function [m, c, pos] = bch1517_decode(r)
%BCH1517_DECODE  Peterson decoder, corrects up to 2 errors.
%   pos = 0-based error positions found ([] if none, NaN if failure).
[EXP, LOG] = gf16_tables();
mul = @(a, b) (a ~= 0 && b ~= 0) * EXP(mod(LOG(a+1) + LOG(b+1), 15) + 1);
dv  = @(a, b) (a ~= 0) * EXP(mod(LOG(a+1) - LOG(b+1), 15) + 1);
pw  = @(a, k) (a ~= 0) * EXP(mod(LOG(a+1) * k, 15) + 1);
s1 = 0; s3 = 0;
for i = 0:14
    if r(i+1)
        s1 = bitxor(s1, EXP(mod(i, 15) + 1));
        s3 = bitxor(s3, EXP(mod(3*i, 15) + 1));
    end
end
c = r; pos = [];
if s1 == 0 && s3 == 0
    m = r(9:15); return
end
if s1 ~= 0 && s3 == pw(s1, 3)
    pos = LOG(s1 + 1);
elseif s1 == 0
    m = r(9:15); pos = NaN; return
else
    sigma2 = dv(bitxor(s3, pw(s1, 3)), s1);
    for j = 0:14                                  % Chien search, x = alpha^-j
        x = EXP(mod(15 - j, 15) + 1);
        if bitxor(bitxor(1, mul(s1, x)), mul(sigma2, mul(x, x))) == 0
            pos(end+1) = j; %#ok<AGROW>
        end
    end
    if numel(pos) ~= 2, m = r(9:15); pos = NaN; return; end
end
c(pos + 1) = 1 - c(pos + 1);
m = c(9:15);
end
