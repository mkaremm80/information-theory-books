"""Demonstrations for Chapter 12 (run: python3 demo.py)."""
import numpy as np
from ldpc import *

rng = np.random.default_rng(12)
# 1. a PEG code and its encoder
code = peg(1008, 504, 3, seed=1)
G, info = generator_from_H(code.H())
print(f"{code.name}: n={code.n} m={code.m} edges={code.E} girth={girth(code)} k={G.shape[0]}")

# 2. one noisy frame: unsatisfied checks after each iteration
fr = np.random.default_rng(1)                          # one frame at Eb/N0 = 2 dB
u = fr.integers(0, 2, (1, 504)); c = encode_G(G, u).astype(np.uint8)
L = awgn_llr(c, 2.0, 0.5, fr)
print("hard-decision bit errors before decoding:", int(((L < 0) != c).sum()),
      " unsatisfied checks:", int(code.syndrome((L < 0).astype(np.uint8)).sum()))
make_layers(code)
for name, f in [("flooding SPA", lambda t: bp_decode(code, L, t, "spa", early_stop=False)[0]),
                ("flooding NMS", lambda t: bp_decode(code, L, t, "nms", early_stop=False)[0]),
                ("layered NMS ", lambda t: layered_decode(code, L, t, "nms")[0])]:
    trace = []
    for t in range(1, 21):
        x = f(t)
        trace.append(int(code.syndrome(x).sum()))
        if trace[-1] == 0:
            break
    print(f"{name}: unsatisfied checks per iteration {trace}  bit errors at end {int((x != c).sum())}")

# 3. a 5G NR code block
nr = NRCode(1024, 2048); make_layers(nr.code, nr.Z)
u = rng.integers(0, 2, (1000, 1024)); cfull, x = nr.encode(u)
for eb in [1.0, 1.5, 2.0]:
    Lr = awgn_llr(x, eb, 0.5, rng)
    xh, it, _ = layered_decode(nr.code, nr.llr_in(Lr), 50, "nms")
    fer = (xh[:, :1024] != u).any(1).mean()
    print(f"{nr.code.name}, Eb/N0 = {eb} dB: FER = {fer:.3f}, average iterations = {it.mean():.1f}")
