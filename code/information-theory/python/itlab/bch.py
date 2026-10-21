"""bch.py -- binary BCH (15,7) code, corrects up to t = 2 errors (optional
advanced extension of Chapter 14).

Field GF(16) built from the primitive polynomial p(x) = x^4 + x + 1, alpha = x.
Generator polynomial g(x) = m1(x) m3(x)
       m1(x) = x^4 + x + 1                 (minimal polynomial of alpha)
       m3(x) = x^4 + x^3 + x^2 + x + 1     (minimal polynomial of alpha^3)
       g(x)  = x^8 + x^7 + x^6 + x^4 + 1
Bit lists hold polynomial coefficients, index i = coefficient of x^i.
Systematic encoding: c(x) = x^8 m(x) + [x^8 m(x) mod g(x)],
so c[0..7] are parity bits and c[8..14] are the 7 message bits.
Decoding (Peterson): syndromes S1 = r(alpha), S3 = r(alpha^3);
  no error  : S1 = 0 (then S3 = 0)
  one error : S3 = S1^3, error at position j with alpha^j = S1
  two errors: error-locator  1 + S1 x + ((S3 + S1^3)/S1) x^2, roots by Chien search.
"""
N, K = 15, 7
G_POLY = [1, 0, 0, 0, 1, 0, 1, 1, 1]          # g(x) = 1 + x^4 + x^6 + x^7 + x^8

# ---------------------------------------------------------------- GF(16)
EXP = [0] * 30
LOG = [0] * 16
v = 1
for i in range(15):
    EXP[i] = v
    LOG[v] = i
    v <<= 1
    if v & 0x10:
        v ^= 0b10011                            # reduce with x^4 + x + 1
for i in range(15, 30):
    EXP[i] = EXP[i - 15]


def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]


def gdiv(a, b):
    if a == 0:
        return 0
    return EXP[(LOG[a] - LOG[b]) % 15]


def gpow(a, k):
    if a == 0:
        return 0
    return EXP[(LOG[a] * k) % 15]


# ---------------------------------------------------------------- encoder
def poly_mod(dividend, divisor):
    r = list(dividend)
    for i in range(len(r) - 1, len(divisor) - 2, -1):
        if r[i]:
            shift = i - (len(divisor) - 1)
            for j, gj in enumerate(divisor):
                r[shift + j] ^= gj
    return r[:len(divisor) - 1]


def encode(m):
    """m: 7 message bits (m[0] = coefficient of x^0). Returns 15 code bits."""
    shifted = [0] * 8 + list(m)                 # x^8 m(x)
    parity = poly_mod(shifted, G_POLY)
    return parity + list(m)


def syndromes(r):
    s1 = s3 = 0
    for i, bit in enumerate(r):
        if bit:
            s1 ^= EXP[i % 15]
            s3 ^= EXP[(3 * i) % 15]
    return s1, s3


def decode(r):
    """Returns (message bits, corrected codeword, error positions found)."""
    r = list(r)
    s1, s3 = syndromes(r)
    if s1 == 0 and s3 == 0:
        return r[8:], r, []
    if s1 != 0 and s3 == gpow(s1, 3):           # single error
        pos = [LOG[s1]]
    elif s1 == 0:                                # >2 errors: cannot correct
        return r[8:], r, None
    else:
        sigma2 = gdiv(s3 ^ gpow(s1, 3), s1)
        pos = []
        for j in range(15):                      # Chien search: test x = alpha^-j
            x = EXP[(15 - j) % 15]
            if (1 ^ gmul(s1, x) ^ gmul(sigma2, gmul(x, x))) == 0:
                pos.append(j)
        if len(pos) != 2:
            return r[8:], r, None                # decoding failure detected
    for j in pos:
        r[j] ^= 1
    return r[8:], r, pos


if __name__ == "__main__":
    m = [1, 0, 1, 1, 0, 0, 1]
    c = encode(m)
    print("message ", m)
    print("codeword", c, " syndromes", syndromes(c))
    r = list(c); r[2] ^= 1; r[11] ^= 1
    mh, ch, pos = decode(r)
    print("received with errors at positions 2 and 11 -> found", pos,
          "message OK:", mh == m)
