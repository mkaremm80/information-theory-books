"""shannon_fano.py -- Shannon-Fano coding (Lab 2).

Conventions
  * symbols are sorted by decreasing probability; ties are broken by symbol;
  * a list is split at the point where |P(left) - P(right)| is smallest
    (the first such point if several are equal);
  * the left (more probable) group gets bit '0', the right group bit '1'.
"""
from math import log2


def shannon_fano(probs):
    """probs: dict symbol -> probability. Returns dict symbol -> codeword."""
    items = sorted(probs.items(), key=lambda kv: (-kv[1], kv[0]))
    codes = {s: "" for s, _ in items}
    steps = []                                  # (group, left, right) for tables

    def split(group):
        if len(group) == 1:
            return
        total = sum(p for _, p in group)
        best_k, best_diff, run = 1, float("inf"), 0.0
        for k in range(1, len(group)):          # left = group[:k]
            run += group[k - 1][1]
            diff = abs(run - (total - run))
            if diff < best_diff - 1e-12:
                best_k, best_diff = k, diff
        left, right = group[:best_k], group[best_k:]
        steps.append(([s for s, _ in left], [s for s, _ in right]))
        for s, _ in left:
            codes[s] += "0"
        for s, _ in right:
            codes[s] += "1"
        split(left)
        split(right)

    if len(items) == 1:                         # one symbol still needs one bit
        return {items[0][0]: "0"}
    split(items)
    shannon_fano.steps = steps
    return codes


def encode(message, codes):
    return "".join(codes[s] for s in message)


def decode(bits, codes):
    """Prefix-free decoding: read bits until they match a codeword."""
    inverse = {c: s for s, c in codes.items()}
    out, cur = [], ""
    for b in bits:
        cur += b
        if cur in inverse:
            out.append(inverse[cur])
            cur = ""
    if cur:
        raise ValueError("bitstream ended inside a codeword")
    return "".join(out)


def average_length(probs, codes):
    """L = sum p(x) * len(code(x)) in bits/symbol."""
    return sum(p * len(codes[s]) for s, p in probs.items())


def report(probs, codes):
    H = -sum(p * log2(p) for p in probs.values() if p > 0)
    L = average_length(probs, codes)
    return {"H": H, "L": L, "efficiency": H / L, "redundancy": 1 - H / L}


if __name__ == "__main__":
    probs = {"A": 0.30, "B": 0.25, "C": 0.15, "D": 0.12, "E": 0.10, "F": 0.08}
    codes = shannon_fano(probs)
    for s in sorted(codes, key=lambda s: -probs[s]):
        print(f"{s}  p = {probs[s]:.2f}  code = {codes[s]}")
    r = report(probs, codes)
    print("H = %.4f  L = %.4f  efficiency = %.2f%%" % (r["H"], r["L"], 100 * r["efficiency"]))
    msg = "ABACADEF"
    bits = encode(msg, codes)
    print("encode", msg, "->", bits, "->", decode(bits, codes))
