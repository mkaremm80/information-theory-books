"""repetition.py -- the (5,1) repetition code (Lab 8).

Each message bit is sent five times. The decoder takes a majority vote of each
5-bit block: d_min = 5, so up to t = floor((5-1)/2) = 2 errors per block are
corrected. Three or more errors in one block cause a decoding error.
Rate R = 1/5: educational and very redundant, not an efficient practical code.
"""
N_REP = 5


def encode(bits, n=N_REP):
    return [b for b in bits for _ in range(n)]


def decode(received, n=N_REP):
    """Majority vote per block. Returns (bits, number of 1s in each block)."""
    out, votes = [], []
    for i in range(0, len(received), n):
        ones = sum(received[i:i + n])
        votes.append(ones)
        out.append(1 if ones > n // 2 else 0)
    return out, votes


def flip(bits, positions):
    """Return a copy of bits with the given 0-based positions inverted."""
    r = list(bits)
    for p in positions:
        r[p] ^= 1
    return r


if __name__ == "__main__":
    msg = [1, 0, 1]
    c = encode(msg)
    print("message", msg, "-> codeword", c)
    r = flip(c, [1, 3, 6, 9])        # 2 errors in block 1, 2 errors in block 2
    print("received", r)
    dec, votes = decode(r)
    print("ones per block", votes, "-> decoded", dec, "OK" if dec == msg else "FAIL")
    r3 = flip(c, [0, 1, 2])          # 3 errors in block 1
    print("3 errors in one block ->", decode(r3)[0], "(failure expected)")
