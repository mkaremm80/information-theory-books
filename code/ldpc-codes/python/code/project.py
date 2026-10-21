"""Final project: a 5G NR LDPC transport-block link with HARQ (Chapter 16).

Transport block -> TB CRC -> base graph selection -> code-block segmentation
(+ CB CRC) -> LDPC encoding -> rate matching (RV) -> BPSK/AWGN -> LLR ->
de-rate-matching (+ HARQ combining) -> layered NMS decoding -> CRC checks.
"""
import json, numpy as np
from ldpc import *

# ---------------------------------------------------------------- CRCs (TS 38.212, 5.1)
CRC_POLY = {"CRC24A": [24, 23, 18, 17, 14, 11, 10, 7, 6, 5, 4, 3, 1, 0],
            "CRC24B": [24, 23, 6, 5, 1, 0],
            "CRC16":  [16, 12, 5, 0]}


class CRC:
    """crc(m) = m P (mod 2). Row i of P is x^(r + k-1-i) mod g(x), built by repeated
    multiplication by x, so the set-up costs only O(k r)."""

    def __init__(self, name, k):
        exps = CRC_POLY[name]
        self.r, self.k = exps[0], k
        g = np.zeros(self.r + 1, np.uint8)
        g[[self.r - e for e in exps]] = 1                 # g[0] = coefficient of x^r
        self.P = np.zeros((k, self.r), np.uint8)
        cur = g[1:].copy()                                 # x^r mod g
        for i in range(k - 1, -1, -1):
            self.P[i] = cur
            top = cur[0]
            cur = np.r_[cur[1:], 0]                        # multiply by x
            if top:
                cur ^= g[1:]

    def attach(self, m):
        return np.concatenate([m, (m.astype(np.int64) @ self.P) % 2], axis=1).astype(np.uint8)

    def check(self, w):
        return ((w[:, :self.k].astype(np.int64) @ self.P) % 2 == w[:, self.k:]).all(axis=1)


def crc_long_division(m, name):
    """Reference CRC by polynomial long division (used to test the class above)."""
    exps = CRC_POLY[name]; r = exps[0]
    g = np.zeros(r + 1, np.uint8); g[[r - e for e in exps]] = 1
    reg = np.r_[np.array(m, np.uint8), np.zeros(r, np.uint8)]
    for i in range(len(m)):
        if reg[i]:
            reg[i:i + r + 1] ^= g
    return reg[-r:]


# ---------------------------------------------------------------- the transport-block link
class TBLink:
    def __init__(self, A, R):
        self.A, self.R = A, R
        self.tbcrc = CRC("CRC24A" if A > 3824 else "CRC16", A)
        B = A + self.tbcrc.r
        self.bg = 2 if (A <= 292 or (A <= 3824 and R <= 0.67) or R <= 0.25) else 1
        Kcb = 8448 if self.bg == 1 else 3840
        self.C = 1 if B <= Kcb else int(np.ceil(B / (Kcb - 24)))
        self.L_cb = 24 if self.C > 1 else 0               # CB CRC only when segmented
        assert B % self.C == 0, "choose A so that B is a multiple of C"
        self.Kp = B // self.C + self.L_cb                  # bits per code block (K')
        self.cbcrc = CRC("CRC24B", self.Kp - 24) if self.C > 1 else None
        self.G = int(round(A / R))                         # coded bits of the transport block
        self.Er = self.G // self.C                         # coded bits per code block
        self.nr = NRCode(self.Kp, 10 ** 6, bg=self.bg)     # mother code, all rows (HARQ)
        make_layers(self.nr.code, self.nr.Z)
        Z = self.nr.Z
        self.N = (66 if self.bg == 1 else 50) * Z          # circular buffer length N_cb = N
        buf = np.arange(2 * Z, 2 * Z + self.N)             # d_k is code bit 2Z + k
        self.buffer = buf[~np.isin(buf, self.nr.filler)]   # filler bits are skipped
        f = (0, 17, 33, 56) if self.bg == 1 else (0, 13, 25, 43)
        self.k0 = {rv: f[rv] * Z for rv in range(4)}       # floor(f N_cb / N) Z with N_cb = N

    def positions(self, rv):
        """Code-bit positions sent with redundancy version rv."""
        start = np.searchsorted(self.buffer, 2 * self.nr.Z + self.k0[rv])
        return self.buffer[(start + np.arange(self.Er)) % len(self.buffer)]

    def encode(self, a):
        b = self.tbcrc.attach(a)
        n_data = self.Kp - self.L_cb
        cws = []
        for r in range(self.C):
            seg = b[:, r * n_data:(r + 1) * n_data]
            cws.append(self.nr.encode(self.cbcrc.attach(seg) if self.cbcrc else seg)[0])
        return cws

    def decode(self, llrs, max_iter=25):
        """llrs: one full-length LLR array per code block. Returns TB-ok flags,
        decoded payload and average iterations per code block."""
        data, ok_cb, its = [], [], 0
        for L in llrs:
            L = L.copy(); L[:, self.nr.filler] = 1e3        # known filler bits
            x, it, _ = layered_decode(self.nr.code, L, max_iter, "nms")
            its = its + it
            cb = x[:, :self.Kp]
            ok_cb.append(self.cbcrc.check(cb) if self.cbcrc else np.ones(len(cb), bool))
            data.append(cb[:, :self.Kp - self.L_cb])
        b = np.concatenate(data, axis=1)
        ok = self.tbcrc.check(b) & np.all(ok_cb, axis=0)
        return ok, b[:, :self.A], its / self.C


def run(A, R, ebn0_list, frames=1000, batch=100, seed=7, rvs=(0, 2)):
    link = TBLink(A, R); rng = np.random.default_rng(seed)
    res = {"ebn0": list(ebn0_list), "tber1": [], "tber2": [], "ber1": [], "it1": []}
    for eb in ebn0_list:
        e1 = e2 = b1 = it1 = n = 0
        while n < frames:
            a = rng.integers(0, 2, (batch, A), dtype=np.uint8)
            cws = link.encode(a)
            soft = [np.zeros((batch, link.nr.n)) for _ in cws]       # HARQ soft buffers
            for t, rv in enumerate(rvs):
                pos = link.positions(rv)
                for i, c in enumerate(cws):
                    np.add.at(soft[i], (slice(None), pos), awgn_llr(c[:, pos], eb, R, rng))
                ok, ahat, its = link.decode(soft)
                if t == 0:
                    e1 += (~ok).sum(); b1 += (ahat != a).sum(); it1 += its.sum()
                else:
                    e2 += (~ok).sum()
            n += batch
        res["tber1"].append(e1 / n); res["tber2"].append(e2 / n)
        res["ber1"].append(b1 / (n * A)); res["it1"].append(it1 / n)
        print(f"A={A:5d} R={R:.3f} Eb/N0={eb:4.2f} dB: TB error 1st tx {e1 / n:.4f}, "
              f"after rv0+rv2 {e2 / n:.4f}, BER {b1 / (n * A):.2e}, avg it {it1 / n:.1f}", flush=True)
    return res


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    for name, k in [("CRC24A", 50), ("CRC24B", 30), ("CRC16", 40)]:
        m = rng.integers(0, 2, (1, k), dtype=np.uint8)
        assert (CRC(name, k).attach(m)[0, k:] == crc_long_division(m[0], name)).all()
    print("CRC24A, CRC24B and CRC16 agree with polynomial long division")
    for A, R in [(12000, 2 / 3), (1000, 1 / 2)]:
        link = TBLink(A, R)
        print(f"A = {A}, R = {R:.3f}: TB CRC {link.tbcrc.r} bits, BG{link.bg}, C = {link.C} code block(s) of "
              f"K' = {link.Kp} bits, Z = {link.nr.Z} (iLS {link.nr.ils}), fillers {link.nr.Kp - link.Kp}, "
              f"G = {link.G}, E_r = {link.Er}, k0 = {list(link.k0.values())}")
        cws = link.encode(rng.integers(0, 2, (2, A), dtype=np.uint8))
        print("   H c^T = 0 for every code block:", all(not link.nr.code.syndrome(c).any() for c in cws),
              "| rv0 and rv2 overlap in", len(np.intersect1d(link.positions(0), link.positions(2))), "positions")
    out = {"A1000": run(1000, 1 / 2, [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0], frames=2000, batch=200),
           "A12000": run(12000, 2 / 3, [1.5, 1.75, 2.0, 2.25, 2.5], frames=400, batch=50),
           "A1000_low": run(1000, 1 / 2, [-3.25, -3.0, -2.75, -2.5, -2.25, -2.0], frames=1000, batch=200)}
    json.dump(out, open("../data/project.json", "w"), indent=1)
