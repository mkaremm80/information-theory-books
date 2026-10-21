"""Monte-Carlo simulations used for the BER/FER figures of the book."""
import json, sys, time, numpy as np
from ldpc import *

def jsonify(r):
    return {k: [float(x) for x in v] if isinstance(v, list) else v for k, v in r.items()}

def run(tag, **kw):
    t = time.time()
    r = simulate(**kw)
    print(tag, f"{time.time() - t:.0f}s", ["%.2e" % f for f in r["fer"]], flush=True)
    r = jsonify(r)
    out[tag.split()[0]] = r
    json.dump(out, open(f"../data/sim_{group}.json", "w"), indent=1)    # save after every curve
    return r

def peg_setup():
    code = peg(1008, 504, 3, seed=1)
    G, info = generator_from_H(code.H())
    make_layers(code)
    enc = lambda u: (lambda c: (c, c))(encode_G(G, u).astype(np.uint8))
    return code, G, info, enc

def nr_setup(K, E):
    nr = NRCode(K, E)
    make_layers(nr.code, nr.Z)
    return nr

E = lambda a, b, s=0.25: list(np.round(np.arange(a, b + 1e-9, s), 2))
group = sys.argv[1]
out = {}
common = dict(max_frames=40000, min_errors=50, seed=2026, stop_fer=2e-4)
if group == "A":          # decoders on PEG(1008, 504)
    code, G, info, enc = peg_setup()
    k, R = 504, 0.5
    for name, rule, lay in [("SPA", "spa", False), ("NMS", "nms", False), ("OMS", "oms", False),
                            ("MS", "ms", False), ("LNMS", "nms", True), ("LSPA", "spa", True)]:
        dec = (lambda L, r=rule: layered_decode(code, L, 50, r)[:2]) if lay else \
              (lambda L, r=rule: bp_decode(code, L, 50, r)[:2])
        hi = 3.0 if rule != "ms" else 3.75
        run(name, encode=enc, decode=dec, k=k, R=R, ebn0_list=E(0.5, hi), info=info,
                        batch=1000, **common)
    dec = lambda L: bit_flip(code, (L < 0).astype(np.uint8), 50)
    run("BF", encode=enc, decode=dec, k=k, R=R, ebn0_list=E(3.0, 7.0, 0.5), info=info,
                    batch=2000, **common)
elif group == "B":        # maximum number of iterations, flooding NMS
    code, G, info, enc = peg_setup()
    for it in [5, 10, 20, 50, 100]:
        dec = lambda L, it=it: bp_decode(code, L, it, "nms")[:2]
        run(f"I{it}", encode=enc, decode=dec, k=504, R=0.5, ebn0_list=E(0.5, 3.0),
                            info=info, batch=1000, **common)
elif group == "C":        # 5G NR, R = 1/2, block lengths
    for K in [128, 256, 512, 1024, 2048, 4096]:
        nr = nr_setup(K, 2 * K)
        dec = lambda L, nr=nr: layered_decode(nr.code, nr.llr_in(L), 50, "nms")[:2]
        hi = {128: 4.0, 256: 3.5, 512: 3.0, 1024: 2.5, 2048: 2.25, 4096: 2.0}[K]
        run(f"K{K} BG{nr.bg} Z{nr.Z}", encode=nr.encode, decode=dec, k=K, R=0.5,
                           ebn0_list=E(0.0, hi), info=np.arange(K), batch=max(200, 100000 // (2 * K)),
                           **(common | ({"max_frames": 20000} if K >= 2048 else {})))
        out[f"K{K}"].update(bg=nr.bg, Z=nr.Z)
elif group == "D":        # 5G NR, K = 1024, code rates
    for R in [1 / 5, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 5 / 6, 8 / 9]:
        Ecw = int(round(1024 / R))
        nr = nr_setup(1024, Ecw)
        dec = lambda L, nr=nr: layered_decode(nr.code, nr.llr_in(L), 50, "nms")[:2]
        lo = {0.2: -1.0, 1/3: -0.5}.get(R, 0.0)
        run(f"R{R:.3f} BG{nr.bg} Z{nr.Z}", encode=nr.encode, decode=dec, k=1024,
                               R=1024 / Ecw, ebn0_list=E(lo, lo + 5.0), info=np.arange(1024),
                               batch=200, **common)
        out[f"R{R:.3f}"].update(bg=nr.bg, Z=nr.Z, E=Ecw)
json.dump(out, open(f"../data/sim_{group}.json", "w"), indent=1)
