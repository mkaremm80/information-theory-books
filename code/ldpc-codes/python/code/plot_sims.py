"""Figures and summary numbers from the Monte-Carlo results (data/sim_*.json)."""
import json, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#008300", "#8a8984", "#1f1f1f"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Inter", "DejaVu Sans"], "mathtext.fontset": "cm",
    "font.size": 10, "axes.edgecolor": "#8a8984", "axes.grid": True, "grid.color": "#e4e3de",
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "savefig.bbox": "tight", "savefig.dpi": 300, "lines.linewidth": 1.8})
OUT, D = "../figs/", "../data/"
LIM = {0.2: -0.964, 1/3: -0.495, 0.5: 0.187, 2/3: 1.059, 0.75: 1.626, 5/6: 2.362, 8/9: 3.033}
summary = open("../out/sim_summary.txt", "w")

def load(g):
    p = D + f"sim_{g}.json"
    return json.load(open(p)) if os.path.exists(p) else {}

def pts(r, key="fer", min_err=10):
    """Points with at least min_err frame errors (reliable estimates)."""
    e, y, f, n = map(np.array, (r["ebn0"], r[key], r["fer"], r.get("frames", [1e9] * len(r["fer"]))))
    ok = f * n >= min_err if "frames" in r else f > 0
    return e[ok], y[ok]

def at(r, target):
    """Eb/N0 where the FER crosses target (log-linear interpolation)."""
    e, f = pts(r, min_err=5)
    for i in range(len(f) - 1):
        if f[i] >= target > f[i + 1]:
            t = (np.log10(target) - np.log10(f[i])) / (np.log10(f[i + 1]) - np.log10(f[i]))
            return e[i] + t * (e[i + 1] - e[i])
    return np.nan

def curve(ax, r, label, color, key="fer", ls="-", marker="o"):
    e, y = pts(r, key)
    ax.semilogy(e, y, ls, marker=marker, ms=3.5, color=color, label=label)

def note(*a):
    print(*a); print(*a, file=summary)

# ---------------------------------------------------------------- A: decoders on PEG(1008,504)
A = load("A")
if A:
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    names = [("SPA", "sum-product (flooding)"), ("LSPA", "sum-product (layered)"), ("NMS", "NMS $\\alpha$=0.75"),
             ("LNMS", "NMS layered"), ("OMS", "OMS $\\beta$=0.5"), ("MS", "min-sum"), ("BF", "bit flipping")]
    for i, (k, lab) in enumerate(names):
        if k in A:
            for j, key in enumerate(["fer", "ber"]):
                curve(ax[j], A[k], lab, S[i], key, ls="--" if k.startswith("L") else "-")
            note(f"A {k}: Eb/N0 at FER 1e-2 = {at(A[k], 1e-2):.2f} dB, at FER 1e-3 = {at(A[k], 1e-3):.2f} dB")
    for j, t in enumerate(["frame error rate", "bit error rate"]):
        ax[j].axvline(0.187, color=S[7], ls=":", lw=1); ax[j].set_xlabel("$E_b/N_0$ (dB)"); ax[j].set_ylabel(t)
        ax[j].set_xlim(0.4, 7.2)
    ax[0].text(0.25, 2e-4, "Shannon\nlimit", fontsize=7, ha="right")
    ax[1].legend(fontsize=7.5, loc="upper right"); ax[0].set_ylim(1e-5, 1.2)
    plt.savefig(OUT + "sim_decoders.pdf"); plt.close()
    fig, a = plt.subplots(figsize=(6, 3.4))
    for i, (k, lab) in enumerate(names[:5]):
        if k in A:
            a.plot(A[k]["ebn0"], A[k]["avg_iter"], "o-" if not k.startswith("L") else "s--", ms=3, color=S[i], label=lab)
    a.set_xlabel("$E_b/N_0$ (dB)"); a.set_ylabel("average iterations"); a.legend(fontsize=8); a.set_ylim(0, 52)
    plt.savefig(OUT + "sim_iters.pdf"); plt.close()
    for k in ["SPA", "NMS", "LNMS", "LSPA"]:
        if k in A:
            note(f"A {k} avg iterations:", dict(zip(A[k]["ebn0"], np.round(A[k]["avg_iter"], 1))))

# ---------------------------------------------------------------- B: iteration limit
B = load("B")
if B:
    fig, a = plt.subplots(figsize=(6, 3.5))
    for i, it in enumerate([5, 10, 20, 50, 100]):
        if f"I{it}" in B:
            curve(a, B[f"I{it}"], f"max {it} iterations", S[i])
            note(f"B I{it}: FER 1e-2 at {at(B[f'I{it}'], 1e-2):.2f} dB, FER 1e-3 at {at(B[f'I{it}'], 1e-3):.2f} dB")
    a.set_xlabel("$E_b/N_0$ (dB)"); a.set_ylabel("frame error rate"); a.legend(fontsize=8); a.set_ylim(1e-4, 1.2)
    plt.savefig(OUT + "sim_maxiter.pdf"); plt.close()

# ---------------------------------------------------------------- C: 5G NR block lengths
C = load("C")
if C:
    fig, a = plt.subplots(figsize=(6.2, 3.7))
    for i, K in enumerate([128, 256, 512, 1024, 2048, 4096]):
        if f"K{K}" in C:
            r = C[f"K{K}"]
            curve(a, r, f"K = {K} (BG{r.get('bg', '?')}, Z = {r.get('Z', '?')})", S[i])
            note(f"C K{K} BG{r.get('bg', '?')} Z{r.get('Z', '?')}: FER 1e-2 at {at(r, 1e-2):.2f} dB, 1e-3 at {at(r, 1e-3):.2f} dB, "
                 f"avg it at last point {r['avg_iter'][-1]:.1f}")
    a.axvline(0.187, color=S[7], ls=":", lw=1); a.text(0.22, 2e-4, "Shannon limit\n(BPSK, R = 1/2)", fontsize=7)
    a.set_xlabel("$E_b/N_0$ (dB)"); a.set_ylabel("frame error rate"); a.legend(fontsize=7.5); a.set_ylim(1e-4, 1.2)
    plt.savefig(OUT + "sim_lengths.pdf"); plt.close()

# ---------------------------------------------------------------- D: 5G NR code rates
Dd = load("D")
if Dd:
    fig, a = plt.subplots(figsize=(6.4, 3.8))
    for i, (k, r) in enumerate(Dd.items()):
        R = float(k[1:]); lab = {0.2: "1/5", 0.333: "1/3", 0.5: "1/2", 0.667: "2/3", 0.75: "3/4", 0.833: "5/6", 0.889: "8/9"}[round(R, 3)]
        curve(a, r, f"R = {lab} (BG{r.get('bg', '?')})", S[i % 8])
        lim = [v for q, v in LIM.items() if abs(q - R) < 0.01][0]
        a.plot([lim], [1.5e-4], "v", color=S[i % 8], ms=6)
        note(f"D R={lab} BG{r.get('bg', '?')} Z{r.get('Z', '?')} E={r.get('E', '?')}: FER 1e-2 at {at(r, 1e-2):.2f} dB, limit {lim:.2f} dB, "
             f"gap {at(r, 1e-2) - lim:.2f} dB")
    a.set_xlabel("$E_b/N_0$ (dB)"); a.set_ylabel("frame error rate"); a.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(1.01, 0.5)); a.set_ylim(1e-4, 1.2)
    a.text(-1.2, 2.2e-4, "triangles: Shannon limits (BPSK)", fontsize=7)
    plt.savefig(OUT + "sim_rates.pdf"); plt.close()

# ---------------------------------------------------------------- E: LDPC versus polar
E = load("E")
if E:
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8), sharey=True)
    for j, N in enumerate([512, 1024]):
        for i, (k, lab) in enumerate([("LDPC_LNMS", "LDPC 5G, layered NMS, 50 it."), ("LDPC_SPA", "LDPC 5G, SPA, 50 it."),
                                      ("POLAR_SC", "polar SC"), ("POLAR_CASCL8", "polar CA-SCL, L = 8, CRC-11")]):
            key = f"N{N}_{k}"
            if key in E:
                curve(ax[j], E[key], lab, S[i], ls="-" if "LDPC" in k else "--")
                note(f"E {key}: FER 1e-2 at {at(E[key], 1e-2):.2f} dB, 1e-3 at {at(E[key], 1e-3):.2f} dB")
        ax[j].set_title(f"N = {N}, K = {N // 2} information bits", fontsize=10); ax[j].set_xlabel("$E_b/N_0$ (dB)")
        ax[j].axvline(0.187, color=S[7], ls=":", lw=1)
    ax[0].set_ylabel("frame error rate"); ax[1].legend(fontsize=8); ax[0].set_ylim(1e-4, 1.2)
    plt.savefig(OUT + "ldpc_vs_polar.pdf"); plt.close()

# ---------------------------------------------------------------- MATLAB versus Python
if os.path.exists("../matlab/sim_matlab.mat") and C:
    import scipy.io as sio
    m = sio.loadmat("../matlab/sim_matlab.mat", simplify_cells=True)["res"]
    fig, a = plt.subplots(figsize=(5.6, 3.4))
    curve(a, C["K1024"], "Python (NumPy)", S[0])
    mm = {k: list(np.atleast_1d(v)) for k, v in m.items()}
    curve(a, mm, "MATLAB (GNU Octave)", S[1], ls="--", marker="s")
    a.set_xlabel("$E_b/N_0$ (dB)"); a.set_ylabel("frame error rate"); a.legend(fontsize=8)
    plt.savefig(OUT + "matlab_vs_python.pdf"); plt.close()

# ---------------------------------------------------------------- final project
if os.path.exists(D + "project.json"):
    P = json.load(open(D + "project.json"))
    low = P.pop("A1000_low", None)
    if low is None and os.path.exists(D + "project_low.json"):
        low = json.load(open(D + "project_low.json"))
    if low:                                              # merge the low-SNR HARQ points
        for k in ["ebn0", "tber1", "tber2", "ber1", "it1"]:
            P["A1000"][k] = list(low[k]) + list(P["A1000"][k])
    fig, a = plt.subplots(figsize=(6.2, 3.6))
    for i, (k, lab) in enumerate([("A1000", "A = 1000, R = 1/2"), ("A12000", "A = 12000, R = 2/3")]):
        r = P[k]; e = np.array(r["ebn0"])
        for key, ls, mk, t in [("tber1", "-", "o", "1st transmission (rv0)"), ("tber2", "--", "s", "after HARQ (rv0 + rv2)")]:
            y = np.array(r[key]); ok = y > 0
            a.semilogy(e[ok], y[ok], ls, marker=mk, ms=3.5, color=S[i], label=f"{lab}: {t}")
        note(f"project {k}: tber1 {r['tber1']}, tber2 {r['tber2']}, it {np.round(r['it1'], 1).tolist()}")
    a.set_xlabel("$E_b/N_0$ of each transmission (dB)"); a.set_ylabel("transport-block error rate")
    a.legend(fontsize=7.5); a.set_ylim(1e-3, 1.2)
    plt.savefig(OUT + "project.pdf"); plt.close()
summary.close()
