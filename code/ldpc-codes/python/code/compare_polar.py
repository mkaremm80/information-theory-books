"""LDPC (5G NR, layered NMS / flooding SPA) versus polar (SC, CA-SCL) at equal
length and rate. polar.py is the library of Part One."""
import json, sys, time, numpy as np
import polar
from ldpc import *

out = {}
def save():
    json.dump(out, open("../data/sim_E.json", "w"), indent=1)

E = lambda a, b: list(np.round(np.arange(a, b + 1e-9, 0.25), 2))
for N in [512, 1024]:
    K = N // 2
    nr = NRCode(K, N); make_layers(nr.code, nr.Z)
    for tag, dec in [("LNMS", lambda L, nr=nr: layered_decode(nr.code, nr.llr_in(L), 50, "nms")[:2]),
                     ("SPA", lambda L, nr=nr: bp_decode(nr.code, nr.llr_in(L), 50, "spa")[:2])]:
        t = time.time()
        r = simulate(nr.encode, dec, K, 0.5, E(0.5, 3.5), info=np.arange(K), max_frames=40000,
                     min_errors=50, batch=500, seed=11, stop_fer=2e-4)
        out[f"N{N}_LDPC_{tag}"] = {k: [float(x) for x in v] for k, v in r.items()} | {"bg": nr.bg, "Z": nr.Z}
        print(f"N{N} LDPC {tag} BG{nr.bg} Z{nr.Z}", f"{time.time()-t:.0f}s", ["%.1e" % f for f in r["fer"]], flush=True)
        save()
    for tag, L, Kp, crc in [("SC", 1, K, None), ("CASCL8", 8, K + 11, "CRC11")]:
        t = time.time(); fer_all = []
        res = {"ebn0": [], "ber": [], "fer": []}
        for eb in E(0.5, 4.0):
            r = polar.simulate(N, Kp, [eb], L=L, crc_name=crc, max_frames=40000, min_errors=50, batch=500, seed=11)
            res["ebn0"].append(eb); res["ber"].append(float(r["ber"][0])); res["fer"].append(float(r["bler"][0]))
            if r["bler"][0] < 2e-4:
                break
        out[f"N{N}_POLAR_{tag}"] = res
        print(f"N{N} polar {tag}", f"{time.time()-t:.0f}s", ["%.1e" % f for f in res["fer"]], flush=True)
        save()
