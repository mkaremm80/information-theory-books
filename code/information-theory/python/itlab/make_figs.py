"""Figures of the Information Theory book (matplotlib part)."""
import random, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import log2
import entropy as E, markov as mk, anomaly as an, channel as chn, compare as cmpm, project as pj

S = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#008300"]
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "DejaVu Sans"],
    "mathtext.fontset": "cm", "font.size": 10, "axes.edgecolor": "#8a8984",
    "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.grid": True,
    "grid.color": "#e4e3de", "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "savefig.bbox": "tight",
    "lines.linewidth": 2})
O = "../figs/"

# 1 binary entropy + min-entropy
p = np.linspace(1e-4, 1 - 1e-4, 500)
fig, ax = plt.subplots(figsize=(6.4, 3.2))
ax.plot(p, [E.binary_entropy(x) for x in p], color=S[0], label=r"Shannon entropy $H_2(p)$")
ax.plot(p, -np.log2(np.maximum(p, 1 - p)), color=S[1], ls="--", label=r"min-entropy $-\log_2\max(p,1-p)$")
for x in (0.11, 0.5):
    ax.plot(x, E.binary_entropy(x), "o", color=S[0], ms=6)
ax.annotate("H(0.11) = 0.50", (0.11, 0.5), xytext=(0.16, 0.25), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=8)
ax.annotate("maximum 1 bit at p = 0.5", (0.5, 1), xytext=(0.55, 1.05), fontsize=8)
ax.set_xlabel("probability p of bit 1"); ax.set_ylabel("bits"); ax.set_ylim(0, 1.15); ax.legend(fontsize=8, loc="lower center")
plt.savefig(O + "binary_entropy.pdf"); plt.close()

# 2 self-information
p = np.linspace(0.005, 1, 400)
fig, ax = plt.subplots(figsize=(6.4, 2.8))
ax.plot(p, -np.log2(p), color=S[3])
for x, lab in [(0.5, "coin: 1 bit"), (1 / 6, "die face: 2.58 bits"), (1 / 64, "1/64: 6 bits")]:
    ax.plot(x, -log2(x), "o", color=S[1]); ax.annotate(lab, (x, -log2(x)), xytext=(8, 4), textcoords="offset points", fontsize=8)
ax.set_xlabel("probability p(x)"); ax.set_ylabel(r"$I(x)=-\log_2 p(x)$ (bits)")
plt.savefig(O + "self_information.pdf"); plt.close()

# 3 password distributions
fig, axs = plt.subplots(1, 2, figsize=(9, 2.8))
uni = [1 / 16] * 16
zipf = np.array([1 / (k + 1) ** 1.2 for k in range(16)]); zipf /= zipf.sum()
for ax, d, t, c in [(axs[0], uni, "uniform: H = 4.00 bits, min-entropy = 4.00 bits", S[0]),
                    (axs[1], zipf, "human-like: H = %.2f bits, min-entropy = %.2f bits" % (E.entropy(list(zipf)), E.min_entropy(list(zipf))), S[1])]:
    ax.bar(range(1, 17), d, color=c); ax.set_title(t, fontsize=9); ax.set_xlabel("password rank"); ax.grid(axis="x", visible=False)
axs[0].set_ylabel("probability")
plt.savefig(O + "password_dists.pdf"); plt.close()
np.save("../data/zipf16.npy", zipf)

# 4 Markov: convergence to stationary + entropy rate vs memoryless
P = [[0.9, 0.1], [0.5, 0.5]]
fig, ax = plt.subplots(figsize=(6.4, 2.8))
pi = [1.0, 0.0]; tr = [pi[0]]
for _ in range(15):
    pi = [pi[0] * P[0][0] + pi[1] * P[1][0], pi[0] * P[0][1] + pi[1] * P[1][1]]; tr.append(pi[0])
pi2 = [0.0, 1.0]; tr2 = [0.0]
for _ in range(15):
    pi2 = [pi2[0] * P[0][0] + pi2[1] * P[1][0], pi2[0] * P[0][1] + pi2[1] * P[1][1]]; tr2.append(pi2[0])
ax.plot(tr, "o-", color=S[0], ms=4, label="start in state 0")
ax.plot(tr2, "s-", color=S[1], ms=4, label="start in state 1")
ax.axhline(5 / 6, color="#8a8984", ls="--", lw=1); ax.text(10, 0.86, r"$\pi_0 = 5/6$", fontsize=9)
ax.set_xlabel("time step n"); ax.set_ylabel(r"$P(X_n = 0)$"); ax.legend(fontsize=8)
plt.savefig(O + "markov_convergence.pdf"); plt.close()

# 5 AEP
rng = np.random.default_rng(3)
p1 = 0.2; Hs = E.binary_entropy(p1)
fig, ax = plt.subplots(figsize=(6.4, 3))
for i, n in enumerate([10, 100, 1000]):
    k = rng.binomial(n, p1, 20000)
    v = -(k * np.log2(p1) + (n - k) * np.log2(1 - p1)) / n
    ax.hist(v, bins=60, density=True, histtype="step", lw=1.8, color=S[i], label=f"n = {n}")
ax.axvline(Hs, color="k", ls="--", lw=1); ax.text(Hs + 0.01, ax.get_ylim()[1] * 0.85, "H = %.3f" % Hs, fontsize=8)
ax.set_xlabel(r"$-\frac{1}{n}\log_2 p(X_1\ldots X_n)$ (bits/symbol)"); ax.set_ylabel("density"); ax.legend(fontsize=8)
plt.savefig(O + "aep.pdf"); plt.close()

# 6 anomaly scores
P_hat, base, res = an.run()
train = mk.simulate(an.P_NORMAL, 20000, seed=7)
vn = an.windows(mk.simulate(an.P_NORMAL, 200 * 50, seed=9), 50)
va = an.windows(mk.simulate(an.P_ATTACK, 200 * 50, seed=10), 50)
sn = [an.kl(an.freq(w), base) for w in vn]; sa = [an.kl(an.freq(w), base) for w in va]
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3))
axs[0].hist(sn, bins=30, color=S[0], alpha=0.8, label="normal windows")
axs[0].hist(sa, bins=30, color=S[1], alpha=0.8, label="attack windows")
axs[0].axvline(res["KL"][0], color="k", ls="--", lw=1.2, label="threshold (99th pct)")
axs[0].set_xlabel("KL score D(window || baseline) (bits)"); axs[0].legend(fontsize=7.5)
th = np.linspace(0, 1.2, 200)
axs[1].plot(th, [np.mean(np.array(sn) > t) for t in th], color=S[0], label="false positive rate")
axs[1].plot(th, [np.mean(np.array(sa) <= t) for t in th], color=S[1], label="false negative rate")
axs[1].set_xlabel("threshold"); axs[1].set_ylabel("rate"); axs[1].legend(fontsize=8)
plt.savefig(O + "anomaly.pdf"); plt.close()

# 7 BSC / BEC capacity
x = np.linspace(0, 1, 300)
fig, ax = plt.subplots(figsize=(6.4, 2.8))
ax.plot(x, [chn.bsc_capacity(v) for v in x], color=S[0], label=r"BSC: $1-H_2(p)$")
ax.plot(x, 1 - x, color=S[2], ls="--", label=r"BEC: $1-\epsilon$")
ax.set_xlabel(r"error / erasure probability ($p$ or $\epsilon$)"); ax.set_ylabel("capacity (bits/use)"); ax.legend(fontsize=8)
plt.savefig(O + "capacity.pdf"); plt.close()

# 8 compression comparison
ds = cmpm.datasets()
names = list(ds); methods = ["Shannon-Fano", "Huffman", "Arithmetic", "LZW"]
R = {n: cmpm.compare(t) for n, t in ds.items()}
fig, ax = plt.subplots(figsize=(9, 3.4))
w = 0.18
for i, m in enumerate(methods):
    vals = [(R[n]["rows"][m][0] + R[n]["rows"][m][1]) / R[n]["n"] for n in names]
    ax.bar(np.arange(len(names)) + (i - 1.5) * w, vals, w, color=S[i], label=m, edgecolor="white")
for j, n in enumerate(names):
    ax.plot([j - 2 * w, j + 2 * w], [R[n]["H"]] * 2, color="k", lw=1.4, ls="--")
ax.set_xticks(range(len(names))); ax.set_xticklabels(names); ax.set_ylabel("bits per symbol (payload + header)")
ax.plot([], [], "k--", lw=1.4, label="zeroth-order entropy H"); ax.legend(fontsize=8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.18))
ax.grid(axis="x", visible=False)
plt.savefig(O + "compare.pdf"); plt.close()
import json; json.dump({n: {"n": R[n]["n"], "H": R[n]["H"], "raw": R[n]["raw"], "rows": R[n]["rows"]} for n in names}, open("../data/compare.json", "w"), indent=1)

# 9 project BER plot
pj.plots(pj.DEFAULT, out_prefix=O + "project")
print("figures ok")
