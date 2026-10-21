"""arithmetic.py -- arithmetic coding (Lab 4).

Two versions are given.

1. encode_exact / decode_exact  (teaching version)
   Uses exact fractions, so every interval in the book can be checked by hand.
   The decoder is told the message length n (or stops at a terminator symbol).

2. encode_int / decode_int  (practical finite-precision version)
   32-bit integer arithmetic with the classic scaling rules
   E1 (both in lower half), E2 (both in upper half) and E3 (middle half,
   'pending' bits) of Witten, Neal and Cleary (1987). Symbol frequencies are
   integers; the decoder is told the message length n.

Conventions: symbols are ordered by sorted symbol name in the cumulative table;
symbol s occupies [cum_low(s), cum_high(s)) inside the current interval.
"""
from fractions import Fraction
from math import ceil, log2


# ------------------------------------------------------------ cumulative table
def cumulative(probs):
    """dict symbol -> (low, high) as exact Fractions, symbols in sorted order."""
    table, c = {}, Fraction(0)
    for s in sorted(probs):
        p = Fraction(probs[s]).limit_denominator(10**9)
        table[s] = (c, c + p)
        c += p
    if c != 1:
        raise ValueError("probabilities must sum to 1 exactly")
    return table


# ------------------------------------------------------------ exact version
def encode_exact(message, probs, trace=False):
    """Returns (low, high, bits): the final interval and the shortest binary
    string whose dyadic interval [0.bits, 0.bits + 2^-len) lies inside it."""
    table = cumulative(probs)
    low, high = Fraction(0), Fraction(1)
    steps = []
    for s in message:
        width = high - low
        c_lo, c_hi = table[s]
        low, high = low + width * c_lo, low + width * c_hi
        steps.append((s, low, high))
    if trace:
        encode_exact.steps = steps
    # shortest tag: try lengths k = 1, 2, ... (k <= ceil(-log2(width)) + 1)
    k = 1
    while True:
        m = ceil(low * 2**k)                     # smallest k-bit number >= low
        if Fraction(m + 1, 2**k) <= high:
            bits = format(m, "0%db" % k)
            return low, high, bits
        k += 1


def decode_exact(bits, probs, n):
    """Decode n symbols from the binary fraction 0.bits."""
    table = cumulative(probs)
    value = Fraction(int(bits, 2), 2**len(bits))
    low, high = Fraction(0), Fraction(1)
    out = []
    for _ in range(n):
        width = high - low
        t = (value - low) / width                # position inside the interval
        for s, (c_lo, c_hi) in table.items():
            if c_lo <= t < c_hi:
                out.append(s)
                low, high = low + width * c_lo, low + width * c_hi
                break
    return "".join(out)


# ------------------------------------------------------------ integer version
PREC = 32
FULL = (1 << PREC) - 1
HALF = 1 << (PREC - 1)
QUARTER = 1 << (PREC - 2)


def _freq_table(counts):
    syms = sorted(counts)
    cum, c = {}, 0
    for s in syms:
        cum[s] = (c, c + counts[s])
        c += counts[s]
    return syms, cum, c


def encode_int(message, counts):
    """counts: dict symbol -> positive integer frequency. Returns a bit string."""
    syms, cum, total = _freq_table(counts)
    if total > QUARTER:
        raise ValueError("total frequency too large for 32-bit precision")
    low, high, pending = 0, FULL, 0
    out = []

    def emit(bit):
        nonlocal pending
        out.append(bit)
        out.extend([1 - bit] * pending)          # resolve pending E3 bits
        pending = 0

    for s in message:
        rng = high - low + 1
        c_lo, c_hi = cum[s]
        high = low + rng * c_hi // total - 1
        low = low + rng * c_lo // total
        while True:
            if high < HALF:                       # E1: interval in lower half
                emit(0)
            elif low >= HALF:                     # E2: interval in upper half
                emit(1)
                low -= HALF; high -= HALF
            elif low >= QUARTER and high < 3 * QUARTER:   # E3: middle half
                pending += 1
                low -= QUARTER; high -= QUARTER
            else:
                break
            low = 2 * low
            high = 2 * high + 1
    pending += 1                                  # flush: two final bits
    emit(0 if low < QUARTER else 1)
    return "".join(map(str, out))


def decode_int(bits, counts, n):
    syms, cum, total = _freq_table(counts)
    stream = [int(b) for b in bits]
    pos = 0

    def next_bit():
        nonlocal pos
        b = stream[pos] if pos < len(stream) else 0   # pad with zeros
        pos += 1
        return b

    low, high, code = 0, FULL, 0
    for _ in range(PREC):
        code = (code << 1) | next_bit()
    out = []
    for _ in range(n):
        rng = high - low + 1
        target = ((code - low + 1) * total - 1) // rng
        for s in syms:                            # find the symbol slot
            c_lo, c_hi = cum[s]
            if c_lo <= target < c_hi:
                break
        out.append(s)
        high = low + rng * c_hi // total - 1
        low = low + rng * c_lo // total
        while True:
            if high < HALF:
                pass
            elif low >= HALF:
                low -= HALF; high -= HALF; code -= HALF
            elif low >= QUARTER and high < 3 * QUARTER:
                low -= QUARTER; high -= QUARTER; code -= QUARTER
            else:
                break
            low = 2 * low
            high = 2 * high + 1
            code = (code << 1) | next_bit()
    return "".join(out)


if __name__ == "__main__":
    probs = {"A": Fraction(1, 2), "B": Fraction(1, 4), "C": Fraction(1, 4)}
    msg = "ABAC"
    low, high, bits = encode_exact(msg, probs, trace=True)
    for s, lo, hi in encode_exact.steps:
        print(f"after '{s}': [{lo}, {hi})  width {hi - lo}")
    print("tag bits:", bits, " -> decoded:", decode_exact(bits, probs, len(msg)))
    text = "TOBEORNOTTOBEORTOBEORNOT"
    from collections import Counter
    counts = Counter(text)
    b = encode_int(text, counts)
    print("integer coder:", len(b), "bits; round trip OK:", decode_int(b, counts, len(text)) == text)
