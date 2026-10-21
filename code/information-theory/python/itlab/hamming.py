"""hamming.py -- Hamming (7,4) and extended Hamming (8,4) SECDED (Lab 7).

Bit positions 1..7 (position-numbered layout of the book):
      position : 1   2   3   4   5   6   7
      content  : p1  p2  d1  p4  d2  d3  d4
  p1 checks positions 1,3,5,7   (positions whose binary form has bit 1 set)
  p2 checks positions 2,3,6,7   (bit 2 set)
  p4 checks positions 4,5,6,7   (bit 4 set)
Even parity is used. The syndrome (s4 s2 s1) read as a binary number is the
position of a single error (0 = no error detected).
"""

def encode74(d):
    """d = [d1, d2, d3, d4] -> 7-bit codeword [p1, p2, d1, p4, d2, d3, d4]."""
    d1, d2, d3, d4 = d
    p1 = d1 ^ d2 ^ d4          # positions 3,5,7
    p2 = d1 ^ d3 ^ d4          # positions 3,6,7
    p4 = d2 ^ d3 ^ d4          # positions 5,6,7
    return [p1, p2, d1, p4, d2, d3, d4]


def syndrome74(r):
    """Returns the error position 0..7 computed from the 7 received bits."""
    s1 = r[0] ^ r[2] ^ r[4] ^ r[6]   # positions 1,3,5,7
    s2 = r[1] ^ r[2] ^ r[5] ^ r[6]   # positions 2,3,6,7
    s4 = r[3] ^ r[4] ^ r[5] ^ r[6]   # positions 4,5,6,7
    return 4 * s4 + 2 * s2 + s1


def decode74(r):
    """Corrects up to one error. Returns (data bits, corrected word, syndrome)."""
    r = list(r)
    s = syndrome74(r)
    if s:
        r[s - 1] ^= 1                 # flip the bit at position s
    return [r[2], r[4], r[5], r[6]], r, s


# Matrix form (same code): codeword c = d G (mod 2), H c^T = 0
G74 = [[1, 1, 1, 0, 0, 0, 0],
       [1, 0, 0, 1, 1, 0, 0],
       [0, 1, 0, 1, 0, 1, 0],
       [1, 1, 0, 1, 0, 0, 1]]
H74 = [[1, 0, 1, 0, 1, 0, 1],     # s1 : column j is the binary form of position j
       [0, 1, 1, 0, 0, 1, 1],     # s2
       [0, 0, 0, 1, 1, 1, 1]]     # s4


# ---------------------------------------------------------------- SECDED (8,4)
def encode84(d):
    """Hamming (7,4) codeword plus an overall parity bit p0 at the end."""
    c = encode74(d)
    return c + [sum(c) % 2]


def decode84(r):
    """Single error correction, double error detection.
    Returns (data or None, status) with status in
    'no error', 'corrected single error', 'parity bit error', 'double error detected'."""
    s = syndrome74(r[:7])
    overall = sum(r) % 2                 # 0 if an even number of bits flipped
    if s == 0 and overall == 0:
        return [r[2], r[4], r[5], r[6]], "no error"
    if overall == 1:                     # odd number of errors: assume one
        if s == 0:
            return [r[2], r[4], r[5], r[6]], "parity bit error"
        d, _, _ = decode74(r[:7])
        return d, "corrected single error"
    return None, "double error detected"  # s != 0 and even parity


if __name__ == "__main__":
    d = [1, 0, 1, 1]
    c = encode74(d)
    print("data", d, "-> codeword", c)
    r = c.copy(); r[4] ^= 1              # error at position 5
    dd, fixed, s = decode74(r)
    print("received", r, "syndrome", s, "-> corrected", fixed, "data", dd)
    r2 = c.copy(); r2[1] ^= 1; r2[5] ^= 1   # two errors (positions 2 and 6)
    dd2, _, s2 = decode74(r2)
    print("two errors -> syndrome", s2, "data", dd2, "(WRONG: Hamming(7,4) corrects only one)")
    c8 = encode84(d)
    r8 = c8.copy(); r8[1] ^= 1; r8[5] ^= 1   # the same two errors in SECDED
    print("SECDED (8,4) codeword", c8, "with two errors ->", decode84(r8)[1])
