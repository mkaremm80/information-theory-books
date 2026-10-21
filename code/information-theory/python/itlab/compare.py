"""compare.py -- fair comparison of the four compressors (Lab 6).

For a text of n symbols:
  raw size            = 8 n bits (8-bit ASCII)
  Shannon-Fano/Huffman = payload bits + header (code table)
  arithmetic (integer) = payload bits + header (frequency table)
  LZW (dynamic widths, 256-entry ASCII start dictionary) = payload bits, no header
Header convention (explicit and simple, used for all table-based coders):
  16 bits for the number of distinct symbols
  + for each symbol: 8 bits (symbol) + 16 bits (its count)
The decoder rebuilds the identical code from the counts, so this header is
sufficient for Shannon-Fano, Huffman and arithmetic coding alike.
"""
from collections import Counter
from math import log2
import shannon_fano as sf
import huffman as hf
import arithmetic as ac
import lzw


def header_bits(counts):
    return 16 + len(counts) * (8 + 16)


def compare(text):
    n = len(text)
    counts = Counter(text)
    probs = {s: c / n for s, c in counts.items()}
    H = -sum(p * log2(p) for p in probs.values())
    rows = {}
    c_sf = sf.shannon_fano(probs)
    b = sf.encode(text, c_sf); assert sf.decode(b, c_sf) == text
    rows["Shannon-Fano"] = (len(b), header_bits(counts))
    c_hf = hf.huffman(probs)
    b = hf.encode(text, c_hf); assert hf.decode(b, c_hf) == text
    rows["Huffman"] = (len(b), header_bits(counts))
    b = ac.encode_int(text, counts); assert ac.decode_int(b, counts, n) == text
    rows["Arithmetic"] = (len(b), header_bits(counts))
    codes = lzw.encode(text); assert lzw.decode(codes) == text
    rows["LZW"] = (lzw.bit_cost(codes, 256), 0)
    return {"n": n, "H": H, "entropy_bits": n * H, "raw": 8 * n, "rows": rows}


DATASETS = {
    "skewed (A-heavy)": "A" * 600 + "B" * 200 + "C" * 120 + "D" * 50 + "E" * 30,
    "English text": None,          # filled in below
    "repetitive": "TOBEORNOT" * 120,
    "uniform random": None,
}

ENGLISH = ("information theory studies how much information a message carries and how "
           "efficiently it can be stored and sent. a source with predictable symbols has "
           "low entropy and can be compressed well. a source with random symbols has high "
           "entropy and cannot be compressed much. compression removes redundancy, and "
           "error correcting codes add redundancy in a controlled way so that a receiver "
           "can find and fix errors. security engineers use the same ideas to measure the "
           "strength of keys and the leakage of side channels. ") * 2


def datasets():
    import random
    rng = random.Random(2026)
    d = dict(DATASETS)
    skew = list(d["skewed (A-heavy)"]); rng.shuffle(skew)
    d["skewed (A-heavy)"] = "".join(skew)
    d["English text"] = ENGLISH
    d["uniform random"] = "".join(rng.choice("ABCDEFGHIJKLMNOP") for _ in range(1000))
    return d


if __name__ == "__main__":
    for name, text in datasets().items():
        r = compare(text)
        print(f"{name}: n = {r['n']}, H = {r['H']:.3f} bits/sym, n*H = {r['entropy_bits']:.0f} bits, raw = {r['raw']}")
        for k, (pay, hdr) in r["rows"].items():
            print(f"   {k:13s} payload {pay:6d}  header {hdr:4d}  total {pay+hdr:6d}  "
                  f"ratio {r['raw']/(pay+hdr):5.2f}  bits/sym {pay/r['n']:.3f}")
