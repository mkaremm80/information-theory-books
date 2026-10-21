"""Numbers for the worked examples of the book (run: python3 examples.py)."""
import itertools, numpy as np
from ldpc import *

np.set_printoptions(precision=4, suppress=True, linewidth=120)
H6 = np.array([[1, 1, 0, 1, 0, 0],
               [0, 1, 1, 0, 1, 0],
               [1, 0, 0, 0, 1, 1],
               [0, 0, 1, 1, 0, 1]], dtype=np.uint8)
code6 = LDPCCode.from_H(H6, "toy(6,3)")

def sec(t):
    print("\n" + "=" * 70 + "\n" + t + "\n" + "=" * 70)

# ------------------------------------------------------------------ chapter 2
sec("CH2 capacity and Shannon limits")
x, w = np.polynomial.hermite_e.hermegauss(120)
def biawgn_cap(sigma):
    y = 1 + sigma * x                       # y given +1
    return 1 - (np.logaddexp(0, -2 * y / sigma ** 2) / np.log(2) * w).sum() / np.sqrt(2 * np.pi)
def limit(R):
    lo, hi = 0.2, 5.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if biawgn_cap(mid) > R else (lo, mid)
    s = lo
    return s, 10 * np.log10(1 / (2 * R * s * s))
for R in [1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 5/6, 8/9]:
    s, eb = limit(R)
    print(f"R={R:.4f}  BI-AWGN sigma*={s:.4f}  Eb/N0 limit={eb:.3f} dB   Gaussian-input limit={10*np.log10((2**R-1)/R):.3f} dB")
print("LLR example y=0.8 sigma^2=0.5 ->", 2 * 0.8 / 0.5)

# ------------------------------------------------------------------ chapter 3
sec("CH3 parity-check matrix")
print("H =\n", H6)
cws = [np.array(c) for c in itertools.product([0, 1], repeat=6) if not ((H6 @ np.array(c)) % 2).any()]
print("codewords:", ["".join(map(str, c)) for c in cws])
print("rank", gf2_rank(H6), "k =", 6 - gf2_rank(H6), "dmin", min(c.sum() for c in cws if c.sum()))
r = np.array([0, 0, 1, 1, 1, 1]); print("r=001111 syndrome", (H6 @ r) % 2)
R_, piv = gf2_rref(H6); print("RREF\n", R_, "pivots", piv)
G, info = generator_from_H(H6); print("G\n", G, "info positions", info)
print("u=101 -> c =", encode_G(G, [1, 0, 1])[0], "check", (H6 @ encode_G(G, [1, 0, 1])[0]) % 2)

# ------------------------------------------------------------------ chapter 4
sec("CH4 Tanner graph and degree distributions")
print("girth toy", girth(code6), "4-cycles", count_4cycles(code6))
lam = {2: 10/30, 3: 12/30, 8: 8/30}; rho = {6: 1.0}
iL = sum(c / d for d, c in lam.items()); iR = sum(c / d for d, c in rho.items())
print("irregular example: int lam", iL, "int rho", iR, "rate", 1 - iR / iL)
print("node perspective L_i:", {d: (c / d) / iL for d, c in lam.items()})
print("lambda(0.5)=", poly(lam, 0.5), "rho(0.5)=", poly(rho, 0.5))
for bg in (1, 2):
    V = nr_base_graph(bg, 0); cd = (V >= 0).sum(0); rd = (V >= 0).sum(1)
    E_ = cd.sum()
    lam5 = {d: (cd == d).sum() * d / E_ for d in np.unique(cd)}
    rho5 = {d: (rd == d).sum() * d / E_ for d in np.unique(rd)}
    print(f"BG{bg}: edges {E_}, VN degree counts", dict(zip(*np.unique(cd, return_counts=True))),
          "CN degree counts", dict(zip(*np.unique(rd, return_counts=True))))
    print("   lambda", {int(k): round(v, 4) for k, v in lam5.items()})

# ------------------------------------------------------------------ chapter 5
sec("CH5 constructions")
for nm, cd in [("Gallager(3,6) n=1008", LDPCCode.from_H(gallager_H(1008, 3, 6, seed=3))),
               ("PEG dv=3 n=1008", peg(1008, 504, 3, seed=1))]:
    print(nm, "rank", gf2_rank(cd.H()), "girth", girth(cd), "4-cycles", count_4cycles(cd))
for Z, base in [(5, [[0, 0, 0, 0], [0, 1, 2, 3]]), (5, [[0, 0, 0, 0], [0, 1, 1, 3]]),
                (7, [[0, 0, 0, 0], [0, 1, 3, 5]])]:
    q = qc_lift(base, Z)
    print("QC base", base, "Z", Z, "girth", girth(q), "4-cycles", count_4cycles(q))
small = peg(12, 6, 2, seed=0); print("PEG(12,6) dv=2 H\n", small.H(), "girth", girth(small))

# ------------------------------------------------------------------ chapter 6
sec("CH6 encoding with a dual-diagonal core")
Zt = 4
base = np.array([[1, 1, -1, 1, 0, -1],
                 [-1, 2, 3, 0, 0, 0],
                 [3, -1, 1, 1, -1, 0]])
toy = qc_lift(base, Zt)
print("toy QC H (12 x 24) rank", gf2_rank(toy.H()), "girth", girth(toy))
s = np.array([1, 0, 0, 1,  0, 1, 1, 0,  1, 0, 1, 0])
H = toy.H().astype(int)
lam_rows = (H[:, :12] @ s) % 2
print("s =", s); print("lambda per block row:", lam_rows.reshape(3, 4))
p0 = lam_rows.reshape(3, 4).sum(0) % 2; print("p0 = sum of lambda blocks =", p0)
# solve remaining blocks by back substitution
def Pm(sh): return np.roll(np.eye(Zt, dtype=int), sh, axis=1)
p1 = (lam_rows[:4] + Pm(1) @ p0) % 2
p2 = (lam_rows[8:] + Pm(1) @ p0) % 2
print("p1 =", p1, " p2 =", p2)
c = np.r_[s, p0, p1, p2]; print("codeword", c, "syndrome zero?", not ((H @ c) % 2).any())
enc = DualEncoder(toy, 12, 12); print("DualEncoder gives same:", (enc.encode(s)[0] == c).all())
cpeg = peg(1008, 504, 3, seed=1); Gp, _ = generator_from_H(cpeg.H())
print("PEG(1008,504): ones in H", cpeg.E, "ones in G", int(Gp.sum()), "density G", Gp.mean())

# ------------------------------------------------------------------ chapter 7
sec("CH7 erasure decoding and bit flipping")
c0 = np.array([0, 0, 1, 0, 1, 1])
y = c0.copy(); y[[0, 2, 5]] = -1
print("BEC received", y)
yy = np.atleast_2d(y).copy()
for rnd in range(1, 5):
    ev = yy[0, code6.vn]
    out, _ = peel(code6, yy)
    print("full peel ->", out[0])
    break
out, rounds = peel(code6, y); print("parallel peeling rounds:", rounds)
yy = np.atleast_2d(y)
for rnd in range(1, 4):
    ev = yy[0, code6.vn]; er = (ev < 0)
    nerr = np.add.reduceat(er.astype(int), code6.cstart)
    print(f" parallel round {rnd}: erasures per check {nerr}")
    yy, _ = peel(code6, yy) if False else (yy, 0)
    break
# step-by-step peeling
yy = y.copy(); step = 0
while (yy < 0).any():
    step += 1; progress = False
    for j in range(4):
        vs = np.nonzero(H6[j])[0]; er = [v for v in vs if yy[v] < 0]
        if len(er) == 1:
            val = sum(yy[v] for v in vs if v != er[0]) % 2
            print(f" round {step}: check {j} (bits {list(vs)}) recovers bit {er[0]} = {val}")
            yy[er[0]] = val; progress = True
    if not progress:
        break
print("result", yy)
# smallest stopping sets
for size in range(1, 7):
    found = []
    for S in itertools.combinations(range(6), size):
        rows = H6[:, list(S)].sum(1)
        if not (rows == 1).any():
            found.append(S)
    if found:
        print("smallest stopping sets size", size, found); break
y2 = c0.copy(); y2[[1, 3]] = -1; print("erase {1,3} ->", peel(code6, y2)[0][0])
from math import erfc, sqrt, log10
Q = lambda x: 0.5 * erfc(x / sqrt(2))
H2 = lambda p: -p * np.log2(p) - (1 - p) * np.log2(1 - p)
for R in [1/2, 1/3, 3/4]:
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = (lo + hi) / 2
        p = Q(sqrt(2 * R * 10 ** (mid / 10)))
        lo, hi = (lo, mid) if 1 - H2(p) > R else (mid, hi)
    print(f"BSC (hard decision) limit R={R:.3f}: Eb/N0 = {hi:.3f} dB, crossover p = {Q(sqrt(2*R*10**(hi/10))):.4f}")
# bit flipping
r = c0.copy(); r[3] ^= 1
print("BF received", r, "syndrome", (H6 @ r) % 2)
fails = H6.T @ ((H6 @ r) % 2); print("failed checks per bit", fails)
x_, it = bit_flip(code6, r, 10); print("BF result", x_[0], "iterations", it)

# ------------------------------------------------------------------ chapter 8
sec("CH8 sum-product example")
yv = np.array([0.9, 1.2, -0.8, -0.3, -1.1, -0.7]); sigma2 = 0.5
L = 2 * yv / sigma2; print("y", yv, "L", L, "hard", (L < 0).astype(int))
v2c = L[code6.vn][None, :]
c2v = check_update(code6, v2c, "spa")
print("edges (check,var):", list(zip(code6.cn, code6.vn)))
print("iteration 1 c2v SPA:", c2v[0])
app = L + code6.var_sum(c2v)[0]; print("APP after it 1:", app, "hard", (app < 0).astype(int))
ch, it, appf = bp_decode(code6, L, 20, "spa"); print("SPA decoded", ch[0], "iterations", it, "APP", appf[0])
ch, it, appf = bp_decode(code6, L, 20, "spa", early_stop=False); print("SPA 20 its APP", appf[0])
# exact MAP
cw = np.array(cws); llik = (np.where(cw == 0, 1, -1) * L / 2).sum(1)  # log p(y|c) up to const
pmap = []
for i in range(6):
    a0 = np.logaddexp.reduce(llik[cw[:, i] == 0]); a1 = np.logaddexp.reduce(llik[cw[:, i] == 1])
    pmap.append(a0 - a1)
print("exact MAP LLRs", np.array(pmap))
Ht = np.array([[1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1]], dtype=np.uint8)
ct = LDPCCode.from_H(Ht, "tree")
cwt = np.array([c for c in itertools.product([0, 1], repeat=7) if not ((Ht @ np.array(c)) % 2).any()])
worst = 0.0
rngt = np.random.default_rng(3)
for trial in range(200):
    Lt = rngt.normal(1.0, 2.0, 7)
    _, _, a_bp = bp_decode(ct, Lt, 10, "spa", early_stop=False)
    ll = (np.where(cwt == 0, 1, -1) * Lt / 2).sum(1)
    mp = np.array([np.logaddexp.reduce(ll[cwt[:, i] == 0]) - np.logaddexp.reduce(ll[cwt[:, i] == 1]) for i in range(7)])
    worst = max(worst, np.abs(a_bp[0] - mp).max())
print("tree code: max |SPA APP - exact MAP| over 200 random LLR vectors =", worst, "girth", girth(ct))
print("v3->c0 message in iteration 2 =", L[3] + c2v[0][10])
print("tanh-rule check: c2v from check 0 to bit 0 =", 2 * np.arctanh(np.tanh(L[1] / 2) * np.tanh(L[3] / 2)))

# ------------------------------------------------------------------ chapter 9
sec("CH9 min-sum variants")
for rule in ["ms", "nms", "oms"]:
    c2v = check_update(code6, L[code6.vn][None, :], rule)
    print(rule, "it1 c2v", c2v[0])
    ch, it, appf = bp_decode(code6, L, 20, rule); print("  decoded", ch[0], "its", it, "APP", appf[0])
make_layers(code6); print("layers of toy code:", [list(code6.cn[e][st]) for e, st, sg in code6.layers])
ch, it, appf = layered_decode(code6, L, 20, "spa"); print("layered SPA", ch[0], it, appf[0])
# iterations comparison on PEG code (same noise)
rng = np.random.default_rng(7); u = rng.integers(0, 2, (2000, 504))
cc = encode_G(Gp, u).astype(np.uint8); make_layers(cpeg)
for eb in [1.5, 2.0, 2.5]:
    Lr = awgn_llr(cc, eb, 0.5, np.random.default_rng(int(eb * 10)))
    r1 = bp_decode(cpeg, Lr, 50, "nms"); r2 = layered_decode(cpeg, Lr, 50, "nms")
    print(f"Eb/N0 {eb}: flooding NMS avg it {r1[1].mean():.2f} FER {(r1[0] != cc).any(1).mean():.4f} | "
          f"layered NMS avg it {r2[1].mean():.2f} FER {(r2[0] != cc).any(1).mean():.4f}  layers {len(cpeg.layers)}")
for a in [0.5, 0.625, 0.7, 0.75, 0.8, 0.875, 1.0]:
    Lr = awgn_llr(cc, 2.0, 0.5, np.random.default_rng(20))
    r1 = bp_decode(cpeg, Lr, 50, "nms", alpha=a)
    print(f"alpha {a}: FER {(r1[0] != cc).any(1).mean():.4f} avg it {r1[1].mean():.2f}")
for b in [0.0, 0.25, 0.5, 0.75, 1.0]:
    Lr = awgn_llr(cc, 2.0, 0.5, np.random.default_rng(20))
    r1 = bp_decode(cpeg, Lr, 50, "oms", beta=b)
    print(f"beta {b}: FER {(r1[0] != cc).any(1).mean():.4f} avg it {r1[1].mean():.2f}")

# ------------------------------------------------------------------ chapter 10
sec("CH10 density evolution")
for eps in [0.40, 0.45]:
    h = de_bec({3: 1.0}, {6: 1.0}, eps, 60)
    print(f"(3,6) eps={eps}: x_l for l=0..10:", np.round(h[:11], 4), "... final", h[-1], "len", len(h))
print("first update by hand eps=.4: 0.4*(1-(1-0.4)^5)^2 =", 0.4 * (1 - 0.6 ** 5) ** 2)
for dv, dc in [(3, 6), (4, 8), (5, 10), (3, 4), (3, 5), (4, 6), (3, 9), (3, 12)]:
    t = threshold_bec({dv: 1.0}, {dc: 1.0}); R = 1 - dv / dc
    print(f"BEC ({dv},{dc}) R={R:.3f} eps*={t:.4f} capacity eps={1-R:.4f} gap={1-R-t:.4f}")
print("stability-type value for dv=2, dc=4: 1/(dc-1) =", 1/3)

# ------------------------------------------------------------------ chapter 11
sec("CH11 5G NR parameters")
for K, E in [(120, 360), (256, 512), (1024, 2048), (1024, 1152), (3840, 19200), (4096, 8192), (8448, 25344)]:
    nr = NRCode(K, E)
    kb = 22 if nr.bg == 1 else (10 if K > 640 else 9 if K > 560 else 8 if K > 192 else 6)
    N = (66 if nr.bg == 1 else 50) * nr.Z
    print(f"K={K:5d} E={E:5d} R={K/E:.3f} BG{nr.bg} kb={kb:2d} Z={nr.Z:3d} iLS={nr.ils} K'={nr.Kp:5d} F={nr.Kp-K:4d} "
          f"N={N:6d} Ncb(after fillers)={nr.Ncb:6d} base rows used={nr.code.m // nr.Z:2d} cols={nr.code.n // nr.Z:2d} edges={nr.code.E}")
print("number of lifting sizes:", sum(len(x) for x in NR_LIFTING_SETS))
for bg in (1, 2):
    V = nr_base_graph(bg, 1)
    print(f"BG{bg}: max code rate with core only = {(22 if bg==1 else 10)}/{(26 if bg==1 else 14)-2} =",
          (22 if bg == 1 else 10) / ((26 if bg == 1 else 14) - 2))
