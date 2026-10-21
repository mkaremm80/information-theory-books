"""Builds the matplotlib figures that do not need the long Monte-Carlo runs."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ldpc import *

OUT = "../figs/"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#008300", "#8a8984"]
INK, GRID = "#1f1f1f", "#e4e3de"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Inter", "DejaVu Sans"],
    "mathtext.fontset": "cm", "font.size": 10, "axes.edgecolor": "#8a8984",
    "axes.labelcolor": INK, "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.dpi": 300, "lines.linewidth": 2})
import sys
PARTS = sys.argv[1:] or ["spy", "limits", "bec", "debec", "exit", "demc"]
log = open("../out/figs_log_" + "_".join(PARTS) + ".txt", "w")
def note(*a):
    print(*a); print(*a, file=log)

# ------------------------------------------------------------------ spy plots
codes = [("Gallager (3,6), n = 1008", LDPCCode.from_H(gallager_H(1008, 3, 6, seed=3))),
         ("PEG, $d_v=3$, n = 1008", peg(1008, 504, 3, seed=1)),
         ("5G NR BG2, Z = 16", qc_lift(np.where(nr_base_graph(2, 0) >= 0, nr_base_graph(2, 0) % 16, -1), 16))]
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
for a, (t, c) in zip(ax, codes):
    a.scatter(c.vn, c.cn, s=0.15, color=S[0] if "NR" not in t else S[1], marker="s", linewidths=0)
    a.set_xlim(0, c.n); a.set_ylim(c.m, 0); a.set_title(t, fontsize=10); a.grid(False)
    a.set_xlabel("column (variable node)"); a.set_aspect("auto")
ax[0].set_ylabel("row (check node)")
plt.savefig(OUT + "spy.pdf"); plt.close()

# ------------------------------------------------------------------ Shannon limits vs rate
x, w = np.polynomial.hermite_e.hermegauss(120)
def cap(sigma):
    y = 1 + sigma * x
    return 1 - (np.logaddexp(0, -2 * y / sigma ** 2) / np.log(2) * w).sum() / np.sqrt(2 * np.pi)
def lim(R):
    lo, hi = 0.1, 8.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if cap(mid) > R else (lo, mid)
    return 10 * np.log10(1 / (2 * R * lo * lo))
Rs = np.linspace(0.05, 0.95, 60)
fig, a = plt.subplots(figsize=(6, 3.4))
a.plot(Rs, [lim(r) for r in Rs], color=S[0], label="BPSK input (BI-AWGN)")
a.plot(Rs, 10 * np.log10((2 ** Rs - 1) / Rs), color=S[1], ls="--", label="Gaussian input (unconstrained)")
for r in [1/5, 1/3, 1/2, 2/3, 5/6]:
    a.plot(r, lim(r), "o", color=S[0], ms=5)
a.set_xlabel("code rate R"); a.set_ylabel("minimum $E_b/N_0$ (dB)"); a.legend()
plt.savefig(OUT + "limits.pdf"); plt.close()

# ------------------------------------------------------------------ BEC peeling simulation
if "bec" in PARTS:
    rng = np.random.default_rng(11)
    eps_grid = np.round(np.arange(0.30, 0.501, 0.01), 2)
    res = {}
    fig, a = plt.subplots(figsize=(6.4, 3.7))
    sets = [("Gallager n = 504", LDPCCode.from_H(gallager_H(504, 3, 6, seed=5)), 2000),
            ("PEG n = 1008", peg(1008, 504, 3, seed=1), 2000),
            ("Gallager n = 2016", LDPCCode.from_H(gallager_H(2016, 3, 6, seed=5)), 2000),
            ("Gallager n = 8064", LDPCCode.from_H(gallager_H(8064, 3, 6, seed=5)), 600)]
    for i, (lab, code, frames) in enumerate(sets):
        H = code.H(); dup = code.n - len(set(map(bytes, H.T)))
        fer = []
        for e in eps_grid:
            bad = 0
            for _ in range(frames // 200):
                y = np.where(rng.random((200, code.n)) < e, -1, 0)     # all-zero codeword (linear code)
                out, _ = peel(code, y)
                bad += (out < 0).any(axis=1).sum()
            fer.append(bad / frames)
        res[lab] = fer
        a.semilogy(eps_grid, np.maximum(fer, 1e-5), "o-", ms=3, color=S[i], label=lab)
        note("BEC peel", lab, "duplicate columns", dup, dict(zip(eps_grid.tolist(), fer)))
    a.plot(eps_grid, eps_grid ** 2, color=S[0], ls=":", lw=1)
    a.text(0.305, 0.05, "$\\epsilon^2$", color=S[0], fontsize=9)
    a.axvline(threshold_bec({3: 1.0}, {6: 1.0}), color=S[3], ls="--", lw=1.4)
    a.axvline(0.5, color=S[6], ls=":", lw=1.4)
    a.text(0.431, 2e-4, "DE threshold\n0.4294", color=S[3], fontsize=8)
    a.text(0.485, 2e-4, "capacity\n0.5", color=S[6], fontsize=8, ha="right")
    a.set_ylim(1e-4, 1.5); a.set_xlabel("erasure probability $\\epsilon$"); a.set_ylabel("frame error rate")
    a.legend(loc="lower left", fontsize=8)
    plt.savefig(OUT + "bec_peel.pdf"); plt.close()
    json.dump(res | {"eps": eps_grid.tolist()}, open("../data/bec_peel.json", "w"))

# ------------------------------------------------------------------ BEC density evolution
lam, rho = {3: 1.0}, {6: 1.0}
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for i, e in enumerate([0.40, 0.42, 0.43, 0.45]):
    h = de_bec(lam, rho, e, 120, tol=1e-9)
    ax[0].semilogy(range(len(h)), np.maximum(h, 1e-9), color=S[i], label=f"$\\epsilon$ = {e}")
ax[0].set_xlabel("iteration $\\ell$"); ax[0].set_ylabel("erasure probability $x_\\ell$"); ax[0].legend()
ax[0].set_ylim(1e-9, 1); ax[0].set_xlim(0, 120)
xs = np.linspace(0, 0.5, 400)
for i, e in enumerate([0.40, 0.4294, 0.45]):
    ax[1].plot(xs, e * poly(lam, 1 - poly(rho, 1 - xs)), color=S[[0, 3, 1][i]], label=f"$f(x)$, $\\epsilon$ = {e}")
ax[1].plot(xs, xs, color=INK, lw=1, ls="--", label="$y = x$")
ax[1].set_xlabel("$x$ (erasure prob. of $v\\to c$ messages)"); ax[1].set_ylabel("$f(x)=\\epsilon\\lambda(1-\\rho(1-x))$")
ax[1].legend(fontsize=8); ax[1].set_xlim(0, 0.5); ax[1].set_ylim(0, 0.5)
plt.savefig(OUT + "de_bec.pdf"); plt.close()

# ------------------------------------------------------------------ EXIT charts
I = np.linspace(0, 1, 400)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
for i, e in enumerate([0.40, 0.45]):
    ax[0].plot(I, 1 - e * (1 - I) ** 2, color=S[i], label=f"VND, $\\epsilon$ = {e}")
ax[0].plot(I ** 5, I, color=INK, lw=1.6, label="CND (inverse), $d_c = 6$")
ax[0].set_xlabel("$I_{A,V}$, $I_{E,C}$"); ax[0].set_ylabel("$I_{E,V}$, $I_{A,C}$"); ax[0].set_title("BEC, (3,6)", fontsize=10)
ax[0].legend(fontsize=8); ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
Jv = np.vectorize(lambda s: J(s)[0]); Jinv = np.vectorize(J_inv)
Ii = np.linspace(1e-4, 1 - 1e-4, 200)
for i, eb in enumerate([1.0, 1.5]):
    R = 0.5; s2ch = 8 * R * 10 ** (eb / 10)
    IEV = Jv(np.sqrt(2 * Jinv(Ii) ** 2 + s2ch))
    ax[1].plot(Ii, IEV, color=S[i], label=f"VND, $E_b/N_0$ = {eb} dB")
IEC = 1 - Jv(np.sqrt(5) * Jinv(1 - Ii))
ax[1].plot(IEC, Ii, color=INK, lw=1.6, label="CND (inverse), $d_c = 6$")
# decoding trajectory at 1.5 dB
s2ch = 8 * 0.5 * 10 ** 0.15; Ia, tx, ty = 0.0, [], []
for _ in range(40):
    Iev = Jv(np.sqrt(2 * Jinv(max(Ia, 1e-6)) ** 2 + s2ch)) if Ia > 0 else Jv(np.sqrt(s2ch))
    tx += [Ia, Ia]; ty += [ty[-1] if ty else 0, float(Iev)]
    Ia = float(1 - Jv(np.sqrt(5) * Jinv(1 - min(Iev, 1 - 1e-9))))
    tx += [Ia]; ty += [float(Iev)]
    if Ia > 0.999:
        break
ax[1].plot(tx, ty, color=S[2], lw=1, label="trajectory, 1.5 dB")
ax[1].set_xlabel("$I_{A,V}$, $I_{E,C}$"); ax[1].set_ylabel("$I_{E,V}$, $I_{A,C}$"); ax[1].set_title("BI-AWGN, (3,6)", fontsize=10)
ax[1].legend(fontsize=8); ax[1].set_xlim(0, 1); ax[1].set_ylim(0, 1)
plt.savefig(OUT + "exit.pdf"); plt.close()
note("EXIT trajectory steps at 1.5 dB:", len(tx) // 3)

# ------------------------------------------------------------------ Monte-Carlo density evolution (AWGN)
fig, a = plt.subplots(figsize=(6.2, 3.5))
for i, sg in enumerate([0.86, 0.875, 0.88, 0.885, 0.90]):
    pe = de_mc(3, 6, sg, iters=300, pop=200000)
    eb = 10 * np.log10(1 / (2 * 0.5 * sg * sg))
    a.semilogy(range(1, len(pe) + 1), np.maximum(pe, 1e-6), color=S[i], label=f"$\\sigma$ = {sg} ({eb:.2f} dB)")
    note("MC DE sigma", sg, "EbN0", round(eb, 3), "iters", len(pe), "final Pe", pe[-1])
a.set_xlabel("iteration"); a.set_ylabel("bit error probability"); a.set_ylim(1e-6, 0.2); a.legend(fontsize=8)
plt.savefig(OUT + "de_mc.pdf"); plt.close()
log.close()
