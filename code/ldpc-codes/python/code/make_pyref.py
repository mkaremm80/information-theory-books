"""Reference results of the Python library, used to test the MATLAB version."""
import numpy as np, scipy.io as sio
from ldpc import *

rng = np.random.default_rng(5)
ref = {}
code = peg(1008, 504, 3, seed=1); H = code.H(); G, info = generator_from_H(H)
ref.update(H_peg=H.astype(float), G_peg=G.astype(float), info_peg=info.astype(float))
u = rng.integers(0, 2, (200, 504)); c = encode_G(G, u).astype(np.uint8)
L = awgn_llr(c, 2.0, 0.5, rng); ref.update(L_peg=L, c_peg=c.astype(float))
make_layers(code)
for rule in ["spa", "ms", "nms", "oms"]:
    x, it, app = bp_decode(code, L, 50, rule)
    ref[f"x_{rule}"] = x; ref[f"it_{rule}"] = it.astype(float)
for rule in ["nms", "spa"]:
    x, it, app = layered_decode(code, L, 50, rule)
    ref[f"xl_{rule}"] = x.astype(float); ref[f"itl_{rule}"] = it.astype(float)
v2c = rng.normal(0, 3, (20, code.E)); ref["v2c"] = v2c
for rule in ["spa", "ms", "nms", "oms"]:
    ref[f"cu_{rule}"] = check_update(code, v2c, rule)
r = (awgn_llr(c[:50], 5.0, 0.5, rng) < 0).astype(np.uint8)
x, it = bit_flip(code, r, 50); ref.update(r_bf=r.astype(float), x_bf=x.astype(float), it_bf=it.astype(float))
y = np.where(rng.random((50, 1008)) < 0.38, -1, c[:50].astype(int))
x, rd = peel(code, y); ref.update(y_peel=y.astype(float), x_peel=x.astype(float), rd_peel=rd.astype(float))
cases = [(40, 200), (120, 360), (500, 1250), (1024, 2048), (1024, 1152), (3840, 19200), (8448, 25344), (100, 3000)]
for i, (K, E) in enumerate(cases):
    nr = NRCode(K, E); uu = rng.integers(0, 2, (3, K)); cc, xx = nr.encode(uu)
    ref[f"nr{i}_KE"] = np.array([K, E], float); ref[f"nr{i}_u"] = uu.astype(float)
    ref[f"nr{i}_c"] = cc.astype(float); ref[f"nr{i}_x"] = xx.astype(float)
    ref[f"nr{i}_par"] = np.array([nr.bg, nr.Z, nr.ils, nr.code.m, nr.code.n, nr.code.E], float)
nr = NRCode(1024, 2048); make_layers(nr.code, nr.Z)
uu = rng.integers(0, 2, (100, 1024)); cc, xx = nr.encode(uu); Lr = awgn_llr(xx, 1.25, 0.5, rng)
xh, it, _ = layered_decode(nr.code, nr.llr_in(Lr), 50, "nms")
ref.update(nr_L=Lr, nr_u=uu.astype(float), nr_xh=xh.astype(float), nr_it=it.astype(float))
ref["th_bec36"] = threshold_bec({3: 1.0}, {6: 1.0}); ref["th_bec34"] = threshold_bec({3: 1.0}, {4: 1.0})
ref["th_ga36"] = threshold_ga({3: 1.0}, {6: 1.0})
for k_, v_ in list(ref.items()):
    if isinstance(v_, np.ndarray) and v_.dtype.kind == "f" and np.isin(v_, [0.0, 1.0]).all() and k_ not in ("it_spa",):
        ref[k_] = v_.astype(np.uint8)                       # binary arrays: compact
    if isinstance(v_, np.ndarray) and v_.dtype.kind == "f" and v_.ndim == 2 and k_ in ("L_peg", "nr_L", "v2c") or k_.startswith("cu_"):
        ref[k_] = ref[k_].astype(np.float64)
sio.savemat("../matlab/pyref_ldpc.mat", ref, do_compression=True)
print("saved", len(ref), "arrays")
