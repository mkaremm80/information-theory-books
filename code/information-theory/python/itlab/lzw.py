"""lzw.py -- Lempel-Ziv-Welch compression (Lab 5).

Conventions
  * the dictionary starts with the single symbols of 'alphabet' with codes
    0, 1, 2, ... in the given order (use alphabet=None for all 256 byte values,
    i.e. code = ASCII value);
  * every output code is followed by one new dictionary entry;
  * with dynamic widths, a code is written with ceil(log2(dictionary size))
    bits, where 'dictionary size' is the size at the moment of writing.
"""
from math import ceil, log2


def initial_dictionary(alphabet=None):
    if alphabet is None:
        return {chr(i): i for i in range(256)}
    return {s: i for i, s in enumerate(alphabet)}


def encode(text, alphabet=None, trace=False):
    """Returns the list of output codes."""
    d = initial_dictionary(alphabet)
    w, out, steps = "", [], []
    for c in text:
        if w + c in d:
            w = w + c                       # keep extending the match
        else:
            out.append(d[w])                # output the longest match
            d[w + c] = len(d)               # add match + next symbol
            steps.append((w, d[w], w + c, len(d) - 1))
            w = c
    if w:
        out.append(d[w])
        steps.append((w, d[w], None, None))
    if trace:
        encode.steps = steps
    return out


def decode(codes, alphabet=None, trace=False):
    """Rebuilds the dictionary while decoding (no dictionary is transmitted)."""
    d = {i: s for s, i in initial_dictionary(alphabet).items()}
    w = d[codes[0]]
    out, steps = [w], [(codes[0], w, None, None)]
    for k in codes[1:]:
        if k in d:
            entry = d[k]
        elif k == len(d):                   # special case: code not yet known
            entry = w + w[0]                # (pattern cScSc)
        else:
            raise ValueError("bad code %d" % k)
        d[len(d)] = w + entry[0]
        steps.append((k, entry, w + entry[0], len(d) - 1))
        out.append(entry)
        w = entry
    if trace:
        decode.steps = steps
    return "".join(out)


def bit_cost(codes, initial_size, dynamic=True, fixed_width=12):
    """Total bits to store the code sequence."""
    if not dynamic:
        return len(codes) * fixed_width
    total, size = 0, initial_size
    for i, _ in enumerate(codes):
        total += max(1, ceil(log2(size)))
        size += 1 if i < len(codes) - 1 else 0   # one new entry per code
    return total


if __name__ == "__main__":
    for text, alpha in [("ABABABA", "AB"), ("TOBEORNOTTOBEORTOBEORNOT", None)]:
        codes = encode(text, alpha, trace=True)
        back = decode(codes, alpha)
        n0 = len(initial_dictionary(alpha))
        print(text)
        print("  codes:", codes)
        print("  round trip OK:", back == text)
        print("  bits (dynamic width):", bit_cost(codes, n0),
              " original:", len(text) * 8, "bits (8-bit ASCII)")
