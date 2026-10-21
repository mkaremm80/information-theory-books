"""project.py -- Final project: an Information Theory & Compression Laboratory.

Pipeline (source coding, then channel coding):
  text -> probabilities -> entropy -> {Shannon-Fano, Huffman, arithmetic, LZW}
       -> comparison table + chart
  Huffman bits -> channel code {none, Hamming(7,4), (5,1) repetition}
       -> binary symmetric channel (or controlled errors) -> channel decoder
       -> Huffman decoder -> recovered text
Run:  python project.py "YOUR TEXT"     (default text if none is given)
"""
import random, sys
from collections import Counter
import entropy as E, shannon_fano as sf, huffman as hf, arithmetic as ac, lzw
import hamming as hm, repetition as rp, channel as chn
from compare import compare


# ---------------------------------------------------------- channel coding layer
def protect(bits, scheme):
    if scheme == "none":
        return list(bits)
    if scheme == "hamming74":
        padded = list(bits) + [0] * (-len(bits) % 4)
        return [b for i in range(0, len(padded), 4) for b in hm.encode74(padded[i:i + 4])]
    if scheme == "rep5":
        return rp.encode(bits)
    raise ValueError(scheme)


def recover(received, scheme, n_bits):
    if scheme == "none":
        return list(received)[:n_bits]
    if scheme == "hamming74":
        out = []
        for i in range(0, len(received), 7):
            out += hm.decode74(received[i:i + 7])[0]
        return out[:n_bits]
    if scheme == "rep5":
        return rp.decode(received)[0][:n_bits]


def safe_huffman_decode(bits, codes):
    """Decode as far as possible (errors may leave a broken tail)."""
    inverse = {c: s for s, c in codes.items()}
    out, cur = [], ""
    for b in bits:
        cur += str(b)
        if cur in inverse:
            out.append(inverse[cur]); cur = ""
    return "".join(out)


def char_errors(a, b):
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def main(text, p=0.01, seed=2026):
    print("=" * 70)
    print("STEP 1-3  source statistics")
    probs = E.probs_from_text(text)
    H = E.entropy(probs)
    print(f"  {len(text)} symbols, {len(probs)} distinct, H = {H:.4f} bits/symbol")
    print("STEP 4-8  source coding and comparison")
    r = compare(text)
    print(f"  {'method':13s}{'payload':>9s}{'header':>8s}{'total':>8s}{'ratio':>7s}")
    for k, (pay, hdr) in r["rows"].items():
        print(f"  {k:13s}{pay:9d}{hdr:8d}{pay + hdr:8d}{r['raw'] / (pay + hdr):7.2f}")
    print(f"  entropy bound n*H = {r['entropy_bits']:.0f} bits, raw = {r['raw']} bits")

    print("STEP 9-12 channel coding over a BSC with p =", p)
    codes = hf.huffman(probs)
    src = [int(b) for b in hf.encode(text, codes)]
    rng = random.Random(seed)
    results = {}
    for scheme in ("none", "hamming74", "rep5"):
        tx = protect(src, scheme)
        rx = [b ^ (rng.random() < p) for b in tx]
        bits = recover(rx, scheme, len(src))
        ber = sum(a != b for a, b in zip(bits, src)) / len(src)
        back = safe_huffman_decode(bits, codes)
        results[scheme] = (len(tx), ber, back == text, char_errors(back, text))
        print(f"  {scheme:10s} sent {len(tx):6d} bits, residual BER {ber:.4f}, "
              f"text recovered exactly: {back == text}, symbol errors {char_errors(back, text)}")

    print("STEP 11   controlled errors: two flips inside one protected block")
    blk = src[:4]
    h = hm.encode74(blk); h[0] ^= 1; h[3] ^= 1
    r5 = rp.encode(blk[:1]); r5[0] ^= 1; r5[3] ^= 1
    print("  Hamming(7,4) with 2 errors ->", hm.decode74(h)[0], "sent", blk,
          "(fails: corrects only 1)")
    print("  (5,1) repetition with 2 errors ->", rp.decode(r5)[0], "sent", blk[:1], "(corrected)")
    return r, results


def plots(text, out_prefix="project"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    r = compare(text)
    names = list(r["rows"])
    tot = [(r["rows"][k][0] + r["rows"][k][1]) / r["n"] for k in names]
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.bar(names, tot, color=["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"])
    ax.axhline(r["H"], color="k", ls="--", lw=1, label="entropy H")
    ax.set_ylabel("bits per symbol (incl. header)"); ax.legend()
    plt.savefig(out_prefix + "_compare.pdf", bbox_inches="tight"); plt.close()
    ps = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
    rng = random.Random(1)
    src = [rng.randint(0, 1) for _ in range(40000)]
    fig, ax = plt.subplots(figsize=(6, 3.2))
    for scheme, lab in [("none", "uncoded"), ("hamming74", "Hamming (7,4)"), ("rep5", "(5,1) repetition")]:
        ber = []
        for p in ps:
            tx = protect(src, scheme)
            rx = [b ^ (rng.random() < p) for b in tx]
            out = recover(rx, scheme, len(src))
            ber.append(max(sum(a != b for a, b in zip(out, src)) / len(src), 1e-6))
        ax.loglog(ps, ber, "o-", label=lab)
    ax.set_xlabel("BSC crossover probability p"); ax.set_ylabel("residual bit error rate")
    ax.legend(); ax.grid(True, which="both", lw=0.3)
    plt.savefig(out_prefix + "_ber.pdf", bbox_inches="tight"); plt.close()


DEFAULT = ("INFORMATION THEORY MEASURES SURPRISE. COMPRESSION REMOVES REDUNDANCY. "
           "ERROR CORRECTION ADDS REDUNDANCY BACK, IN A CONTROLLED WAY. ") * 4

if __name__ == "__main__":
    text = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    main(text)
