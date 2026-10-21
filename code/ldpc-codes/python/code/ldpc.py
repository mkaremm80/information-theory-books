"""ldpc.py -- a small, complete LDPC library for the book
"LDPC Codes: From Zero to Expert" (Part Two).  Only NumPy is needed.

Conventions (same as Part One): bits are 0/1 integers, BPSK maps 0 -> +1 and
1 -> -1, and an LLR is L = ln P(bit=0|y) / P(bit=1|y), so L > 0 favours 0.
"""
import numpy as np
from pathlib import Path

# ==== [gf2]
# ---------------------------------------------------------------- 1. GF(2) algebra
def gf2_rref(A):
    """Reduced row echelon form over GF(2). Returns (R, pivot_columns)."""
    R = np.array(A, dtype=np.uint8) % 2
    m, n = R.shape
    pivots, row = [], 0
    for col in range(n):
        if row == m:
            break
        hits = np.nonzero(R[row:, col])[0]
        if hits.size == 0:
            continue
        p = row + hits[0]
        R[[row, p]] = R[[p, row]]                   # swap the pivot row up
        others = np.nonzero(R[:, col])[0]
        others = others[others != row]
        R[others] ^= R[row]                         # clear the column
        pivots.append(col)
        row += 1
    return R[:row], pivots


def gf2_rank(A):
    return len(gf2_rref(A)[1])


def gf2_inv(A):
    """Inverse of a square, invertible GF(2) matrix."""
    n = A.shape[0]
    R, piv = gf2_rref(np.hstack([A % 2, np.eye(n, dtype=np.uint8)]))
    if piv[:n] != list(range(n)):
        raise ValueError("matrix is singular over GF(2)")
    return R[:, n:]


def generator_from_H(H):
    """Systematic-style generator matrix G (k x n) with H G^T = 0.
    Returns (G, info_pos): codeword bits at info_pos are the message bits."""
    R, piv = gf2_rref(H)
    n = H.shape[1]
    free = [j for j in range(n) if j not in set(piv)]
    G = np.zeros((len(free), n), dtype=np.uint8)
    G[:, free] = np.eye(len(free), dtype=np.uint8)
    G[:, piv] = R[:, free].T                        # c[piv_i] = sum_j R[i, free_j] c[free_j]
    return G, np.array(free)
# ==== [/gf2]


# ==== [code]
# ---------------------------------------------------------------- 2. the code object
class LDPCCode:
    """An LDPC code stored as the edge list of its Tanner graph.
    Edge e joins check node cn[e] and variable node vn[e]; edges are sorted
    by check node (then variable node)."""

    def __init__(self, m, n, cn, vn, name="LDPC"):
        order = np.lexsort((vn, cn))
        self.m, self.n, self.name = m, n, name
        self.cn, self.vn = np.asarray(cn)[order], np.asarray(vn)[order]
        self.E = len(self.cn)
        self.dc = np.bincount(self.cn, minlength=m)          # check degrees
        self.dv = np.bincount(self.vn, minlength=n)          # variable degrees
        self.cstart = np.r_[0, np.cumsum(self.dc)[:-1]]     # first edge of each check
        self.vperm = np.lexsort((self.cn, self.vn))         # edges sorted by variable
        self.vstart = np.r_[0, np.cumsum(self.dv)[:-1]]
        self.layers = None

    @classmethod
    def from_H(cls, H, name="LDPC"):
        cn, vn = np.nonzero(np.asarray(H) % 2)
        return cls(H.shape[0], H.shape[1], cn, vn, name)

    def H(self):
        """Dense parity-check matrix (use only for small codes)."""
        H = np.zeros((self.m, self.n), dtype=np.uint8)
        H[self.cn, self.vn] = 1
        return H

    def design_rate(self):
        return 1 - self.m / self.n

    def check_sum(self, x):
        """Per-check sums of edge values x (shape (B, E)) -> (B, m)."""
        return np.add.reduceat(x, self.cstart, axis=1)

    def var_sum(self, x):
        """Per-variable sums of edge values x (shape (B, E)) -> (B, n)."""
        return np.add.reduceat(x[:, self.vperm], self.vstart, axis=1)

    def syndrome(self, c):
        """Syndrome s = H c^T for a batch of words c (shape (B, n))."""
        c = np.atleast_2d(c)
        return self.check_sum(c[:, self.vn].astype(np.int32)) % 2
# ==== [/code]


# ==== [construct]
# ---------------------------------------------------------------- 3. constructions
def gallager_H(n, dv, dc, seed=0):
    """Gallager's (dv, dc)-regular construction: dv stacked bands, each a
    random column permutation of the first band."""
    rng = np.random.default_rng(seed)
    mb = n // dc                                   # checks per band
    band = np.kron(np.eye(mb, dtype=np.uint8), np.ones((1, dc), dtype=np.uint8))
    bands = [band] + [band[:, rng.permutation(n)] for _ in range(dv - 1)]
    return np.vstack(bands)


def peg(n, m, vn_degrees, seed=0):
    """Progressive edge-growth (Hu, Eleftheriou and Arnold, 2005).
    Each new edge goes to a check node as far as possible from the variable."""
    rng = np.random.default_rng(seed)
    vdeg = np.broadcast_to(vn_degrees, (n,))
    chk_nb = [[] for _ in range(m)]                 # neighbours of each check
    var_nb = [[] for _ in range(n)]
    cdeg = np.zeros(m, dtype=int)
    for j in np.argsort(vdeg, kind="stable"):      # low degrees first
        for k in range(vdeg[j]):
            if k == 0:
                cand = np.flatnonzero(cdeg == cdeg.min())
            else:                                   # breadth-first search from j
                seen_c = np.zeros(m, bool); seen_c[var_nb[j]] = True
                frontier, cand = list(var_nb[j]), None
                while True:
                    nxt_v = {v for c in frontier for v in chk_nb[c]}
                    nxt_c = {c for v in nxt_v for c in var_nb[v]} - set(np.flatnonzero(seen_c))
                    if not nxt_c or seen_c.all():  # cannot grow: take the unreached
                        cand = np.flatnonzero(~seen_c)
                        break
                    seen_c[list(nxt_c)] = True
                    frontier = list(nxt_c)
                    if seen_c.all():               # all reached: farthest layer
                        cand = np.array(sorted(nxt_c))
                        break
                cand = cand[cdeg[cand] == cdeg[cand].min()]
            c = rng.choice(cand)
            chk_nb[c].append(j); var_nb[j].append(c); cdeg[c] += 1
    cn = np.array([c for j in range(n) for c in var_nb[j]])
    vn = np.array([j for j in range(n) for _ in var_nb[j]])
    return LDPCCode(m, n, cn, vn, f"PEG({n},{n - m})")


# ==== [qclift]
def qc_lift(base, Z):
    """Expand a base (exponent) matrix: -1 -> Z x Z zero block, s >= 0 ->
    identity cyclically shifted to the right by s (row r has a 1 in column (r+s) mod Z)."""
    base = np.asarray(base)
    r = np.arange(Z)
    cn, vn = [], []
    for i, j in zip(*np.nonzero(base >= 0)):
        cn.append(i * Z + r)
        vn.append(j * Z + (r + base[i, j]) % Z)
    return LDPCCode(base.shape[0] * Z, base.shape[1] * Z,
                    np.concatenate(cn), np.concatenate(vn), f"QC(Z={Z})")
# ==== [/qclift]


def girth(code):
    """Length of the shortest cycle of the Tanner graph (BFS from each variable)."""
    vnb = [code.cn[code.vperm[s:s + d]] for s, d in zip(code.vstart, code.dv)]
    cnb = [code.vn[s:s + d] for s, d in zip(code.cstart, code.dc)]
    best = np.inf
    for root in range(code.n):
        dist = {("v", root): 0}; parent = {("v", root): None}
        queue = [("v", root)]
        for node in queue:
            kind, idx = node
            nbrs = [("c", c) for c in vnb[idx]] if kind == "v" else [("v", v) for v in cnb[idx]]
            for nb in nbrs:
                if nb == parent[node]:
                    continue
                if nb in dist:                       # a cycle closes here
                    best = min(best, dist[node] + dist[nb] + 1)
                else:
                    dist[nb] = dist[node] + 1; parent[nb] = node; queue.append(nb)
            if 2 * dist[node] >= best:
                break
    return best


def count_4cycles(code):
    """Number of 4-cycles = pairs of checks sharing two (or more) variables."""
    H = code.H().astype(np.int64)
    overlap = H @ H.T
    iu = np.triu_indices(code.m, 1)
    k = overlap[iu]
    return int(np.sum(k * (k - 1) // 2))
# ==== [/construct]


# ==== [nr5g]
# ---------------------------------------------------------------- 4. 5G NR base graphs
NR_LIFTING_SETS = [[2, 4, 8, 16, 32, 64, 128, 256], [3, 6, 12, 24, 48, 96, 192, 384],
                   [5, 10, 20, 40, 80, 160, 320], [7, 14, 28, 56, 112, 224],
                   [9, 18, 36, 72, 144, 288], [11, 22, 44, 88, 176, 352],
                   [13, 26, 52, 104, 208], [15, 30, 60, 120, 240]]   # TS 38.212 Table 5.3.2-1
DATA = Path(__file__).resolve().parent.parent / "data"


def nr_base_graph(bg, ils):
    """Base matrix V (entries -1 or V_ij) of BG1 (46 x 68) or BG2 (42 x 52) for
    lifting-set index ils, read from TS 38.212 Tables 5.3.2-2/3 (csv)."""
    shape = (46, 68) if bg == 1 else (42, 52)
    V = -np.ones(shape, dtype=int)
    row = 0
    for line in open(DATA / f"5G_bg{bg}.csv").read().splitlines()[2:]:
        f = line.split(";")
        if f[0].strip():
            row = int(f[0])
        V[row, int(f[1])] = int(f[2 + ils])
    return V


def nr_select(K, R, bg=None):
    """Base graph and lifting size for K information bits (incl. CRC) and rate R
    (TS 38.212, Sections 7.2.2 and 5.2.2)."""
    if bg is None:
        bg = 2 if (K <= 292 or (K <= 3824 and R <= 0.67) or R <= 0.25) else 1
    if bg == 1:
        kb = 22
    else:
        kb = 10 if K > 640 else 9 if K > 560 else 8 if K > 192 else 6
    Z = min(z for s in NR_LIFTING_SETS for z in s if kb * z >= K)
    ils = next(i for i, s in enumerate(NR_LIFTING_SETS) if Z in s)
    return bg, Z, ils


class NRCode:
    """A 5G NR LDPC code: lifted base graph, encoder and simple rate matching
    (redundancy version 0, no interleaving). Only the rows of the base graph that
    are needed for the E transmitted bits are kept (the other parity bits are not
    sent, so their checks carry no information)."""

    def __init__(self, K, E, bg=None):
        self.bg, Z, ils = nr_select(K, K / E, bg)
        self.K, self.E, self.Z, self.ils = K, E, Z, ils
        self.V = nr_base_graph(self.bg, ils)
        self.P = np.where(self.V >= 0, self.V % Z, -1)      # shift values P = V mod Z
        self.kb = 22 if self.bg == 1 else 10                 # information columns
        self.Kp = self.kb * Z                                # K' incl. filler bits
        self.filler = np.arange(K, self.Kp)                  # known zero bits
        nfull = self.V.shape[1] * Z                          # 68Z or 52Z
        tx = np.setdiff1d(np.arange(2 * Z, nfull), self.filler)   # first 2Z punctured
        self.Ncb = len(tx)
        tx = tx[:E] if E <= len(tx) else tx                   # bits that are ever sent
        ncol = max(self.kb + 4, -(-(tx.max() + 1) // Z))      # base columns in use
        self.code = qc_lift(self.P[:ncol - self.kb, :ncol], Z)
        self.code.name = f"NR BG{self.bg} K={K} E={E} Z={Z}"
        self.n = ncol * Z
        self.enc = DualEncoder(self.code, self.Kp, 4 * Z)
        self.tx_pos0 = tx
        self.tx_pos = np.resize(tx, E)                       # circular buffer, rv = 0

    def encode(self, u):
        u = np.atleast_2d(u)
        s = np.zeros((u.shape[0], self.Kp), dtype=np.uint8)
        s[:, :self.K] = u                                    # filler bits stay 0
        c = self.enc.encode(s)
        return c, c[:, self.tx_pos]                          # full codeword, sent bits

    def llr_in(self, llr_tx):
        """Map received LLRs of the E sent bits to the full-length decoder input."""
        L = np.zeros((llr_tx.shape[0], self.n))
        for p in range(0, self.E, len(self.tx_pos0)):         # repeated bits add up
            q = min(p + len(self.tx_pos0), self.E)
            L[:, self.tx_pos[p:q]] += llr_tx[:, p:q]
        L[:, self.filler] = 1e3                              # filler bits are known zeros
        return L
# ==== [/nr5g]


# ==== [encode]
# ---------------------------------------------------------------- 5. encoders
def encode_G(G, u):
    """Encoding with a generator matrix: c = u G (mod 2)."""
    return (np.atleast_2d(u).astype(np.int64) @ G) % 2


class DualEncoder:
    """Linear-time-style encoder for H = [[A, B, 0], [C1, C2, I]]:
    columns = k systematic bits | g core parity bits | (m-g) extension parity bits.
    Step 1: lam = A s.  Step 2: p1 = B^-1 lam.  Step 3: p2 = C1 s + C2 p1."""

    def __init__(self, code, k, g):
        self.code, self.k, self.g = code, k, g
        H = np.zeros((g, k + g), dtype=np.uint8)
        sel = (code.cn < g) & (code.vn < k + g)
        H[code.cn[sel], code.vn[sel]] = 1
        self.Binv = gf2_inv(H[:, k:k + g])                  # small g x g inverse
        diag = (code.cn >= g) & (code.vn == code.cn + k)
        self.ext_edges = (code.cn >= g) & ~diag             # edges of [C1 C2]
        self.core_edges = (code.cn < g) & (code.vn < k)     # edges of A

    def _rows(self, mask, c, rows):
        """Parity of the bits c[vn] over the selected edges, per check row."""
        B, m = c.shape[0], self.code.m
        idx = (np.arange(B)[:, None] * m + self.code.cn[mask][None, :]).ravel()
        out = np.bincount(idx, weights=c[:, self.code.vn[mask]].ravel(), minlength=B * m)
        return out.reshape(B, m)[:, rows].astype(np.int64) % 2

    def encode(self, s):
        s = np.atleast_2d(s).astype(np.int64)
        B, n = s.shape[0], self.code.n
        c = np.zeros((B, n), dtype=np.int64)
        c[:, :self.k] = s
        lam = self._rows(self.core_edges, c, np.arange(self.g))         # step 1
        c[:, self.k:self.k + self.g] = (lam @ self.Binv.T.astype(np.int64)) % 2   # step 2
        c[:, self.k + self.g:] = self._rows(self.ext_edges, c, np.arange(self.g, self.code.m))  # step 3
        return c.astype(np.uint8)
# ==== [/encode]


# ==== [channel]
# ---------------------------------------------------------------- 6. channels
def bpsk(c):
    return 1.0 - 2.0 * c


def sigma_from_ebn0(ebn0_db, R):
    """Noise standard deviation for BPSK (Es = 1) at the given Eb/N0 and rate."""
    return np.sqrt(1.0 / (2 * R * 10 ** (ebn0_db / 10)))


def awgn_llr(c, ebn0_db, R, rng):
    """Send codewords c through BPSK + AWGN and return channel LLRs 2y/sigma^2."""
    sigma = sigma_from_ebn0(ebn0_db, R)
    y = bpsk(c) + sigma * rng.standard_normal(c.shape)
    return 2 * y / sigma ** 2
# ==== [/channel]


# ==== [hard]
# ---------------------------------------------------------------- 7. hard-decision decoders
def bit_flip(code, r, max_iter=50):
    """Gallager-style parallel bit flipping on hard decisions r (B, n): in each
    iteration flip every bit that sits in the largest number of failed checks."""
    x = np.atleast_2d(r).astype(np.uint8).copy()
    iters = np.zeros(x.shape[0], dtype=int)
    for it in range(1, max_iter + 1):
        s = code.syndrome(x)                                  # (B, m)
        active = s.any(axis=1)
        if not active.any():
            break
        iters[active] = it
        fails = code.var_sum(s[:, code.cn])                   # failed checks per bit
        worst = fails.max(axis=1, keepdims=True)
        flip = (fails == worst) & (worst > 0) & active[:, None]
        x ^= flip.astype(np.uint8)
    return x, iters


def peel(code, y):
    """Peeling decoder for the binary erasure channel. y holds 0, 1 or -1 (erasure).
    A check with exactly one erased neighbour recovers it. Returns the word
    (remaining erasures = -1) and the number of rounds used."""
    y = np.atleast_2d(y).astype(np.int64).copy()
    rounds = np.zeros(y.shape[0], dtype=int)
    for it in range(1, code.n + 1):
        ev = y[:, code.vn]                                    # edge values
        n_er = code.check_sum((ev < 0).astype(np.int64))      # erasures per check
        par = code.check_sum(np.where(ev < 0, 0, ev)) % 2     # parity of known bits
        solve = n_er == 1
        if not solve.any():
            break
        e_solve = solve[:, code.cn] & (ev < 0)                # the single erased edge
        b, e = np.nonzero(e_solve)
        y[b, code.vn[e]] = par[b, code.cn[e]]
        rounds[np.unique(b)] = it
    return y, rounds
# ==== [/hard]


# ==== [bp]
# ---------------------------------------------------------------- 8. soft decoders (flooding)
def phi(x):
    """Gallager's function phi(x) = -ln tanh(x/2); it is its own inverse."""
    x = np.clip(x, 1e-12, 50.0)
    return -np.log(np.tanh(x / 2))


def check_update(code, v2c, rule="spa", alpha=0.75, beta=0.5, starts=None, seg=None):
    """Check-node messages from variable-to-check messages v2c (B, E).
    rule: 'spa' (sum-product), 'ms' (min-sum), 'nms' (normalised, x alpha),
    'oms' (offset, - beta). starts/seg describe the checks (default: all)."""
    st = code.cstart if starts is None else starts
    seg = code.cn if seg is None else seg                         # check of each edge
    neg = (v2c < 0).astype(np.int64)
    sign_all = np.add.reduceat(neg, st, axis=1) % 2              # parity of negatives
    sign = 1 - 2 * ((sign_all[:, seg] + neg) % 2)                 # exclude own sign
    mag = np.abs(v2c)
    if rule == "spa":
        p = phi(mag)
        out = phi(np.add.reduceat(p, st, axis=1)[:, seg] - p)
    else:
        m1 = np.minimum.reduceat(mag, st, axis=1)                  # smallest magnitude
        is1 = mag == m1[:, seg]
        n1 = np.add.reduceat(is1.astype(np.int64), st, axis=1)     # how many equal it
        m2 = np.minimum.reduceat(np.where(is1, np.inf, mag), st, axis=1)
        m2 = np.where(n1 >= 2, m1, m2)                             # second smallest
        out = np.where(is1, m2[:, seg], m1[:, seg])
        if rule == "nms":
            out = alpha * out
        elif rule == "oms":
            out = np.maximum(out - beta, 0.0)
    return sign * out


def bp_decode(code, llr, max_iter=50, rule="spa", alpha=0.75, beta=0.5, early_stop=True):
    """Flooding belief propagation. Returns hard decisions (B, n), the number of
    iterations used per frame and the final a-posteriori LLRs."""
    llr = np.atleast_2d(llr).astype(float)
    B = llr.shape[0]
    app, iters = llr.copy(), np.full(B, max_iter)
    act = np.arange(B)                                    # frames still decoding
    c2v = np.zeros((B, code.E))
    for it in range(1, max_iter + 1):
        v2c = app[act][:, code.vn] - c2v                  # extrinsic: remove own message
        c2v = check_update(code, v2c, rule, alpha, beta)
        app[act] = llr[act] + code.var_sum(c2v)           # a-posteriori LLRs
        if early_stop:
            ok = ~code.syndrome((app[act] < 0).astype(np.uint8)).any(axis=1)
            iters[act[ok]] = it
            act, c2v = act[~ok], c2v[~ok]                 # stop decoding valid codewords
            if act.size == 0:
                break
    return (app < 0).astype(np.uint8), iters, app
# ==== [/bp]


# ==== [layered]
# ---------------------------------------------------------------- 9. layered decoding
def make_layers(code, Z=None):
    """Split the checks into layers whose checks share no variable. For a QC code
    pass Z: every block row of Z checks is one layer."""
    if Z is not None:
        groups = [np.arange(r, r + Z) for r in range(0, code.m, Z)]
    else:                                                # greedy colouring of checks
        used, groups = [], []
        for c in range(code.m):
            vs = set(code.vn[code.cstart[c]:code.cstart[c] + code.dc[c]])
            for g, u in zip(groups, used):
                if not (vs & u):
                    g.append(c); u |= vs; break
            else:
                groups.append([c]); used.append(set(vs))
        groups = [np.array(g) for g in groups]
    layers = []
    for g in groups:
        e = np.concatenate([np.arange(code.cstart[c], code.cstart[c] + code.dc[c]) for c in g])
        layers.append((e, np.r_[0, np.cumsum(code.dc[g])[:-1]], np.repeat(np.arange(len(g)), code.dc[g])))
    code.layers = layers
    return layers


def layered_decode(code, llr, max_iter=50, rule="nms", alpha=0.75, beta=0.5):
    """Row-layered (turbo-decoding message passing) decoder: the a-posteriori
    LLRs are updated after every layer, so new information is used at once."""
    if code.layers is None:
        make_layers(code)
    app = np.atleast_2d(llr).astype(float).copy()
    B = app.shape[0]
    iters = np.full(B, max_iter)
    act = np.arange(B)
    a, c2v = app.copy(), np.zeros((B, code.E))          # working copies of active frames
    for it in range(1, max_iter + 1):
        for e, st, sg in code.layers:
            v2c = a[:, code.vn[e]] - c2v[:, e]              # remove the old message
            new = check_update(code, v2c, rule, alpha, beta, starts=st, seg=sg)
            a[:, code.vn[e]] = v2c + new                    # update APP immediately
            c2v[:, e] = new
        ok = ~code.syndrome((a < 0).astype(np.uint8)).any(axis=1)
        app[act] = a
        iters[act[ok]] = it
        act, a, c2v = act[~ok], a[~ok], c2v[~ok]
        if act.size == 0:
            break
    return (app < 0).astype(np.uint8), iters, app
# ==== [/layered]


# ==== [de]
# ---------------------------------------------------------------- 10. density evolution and EXIT
def poly(coef, x):
    """Evaluate a degree distribution {degree: fraction of edges} as sum c x^(d-1)."""
    return sum(c * x ** (d - 1) for d, c in coef.items())


def de_bec(lam, rho, eps, iters=2000, tol=1e-12):
    """Density evolution on the BEC: x_{l+1} = eps * lam(1 - rho(1 - x_l)).
    Returns the list of erasure probabilities x_l."""
    x, hist = eps, [eps]
    for _ in range(iters):
        x = eps * poly(lam, 1 - poly(rho, 1 - x))
        hist.append(x)
        if x < tol:
            break
    return hist


def threshold_bec(lam, rho, iters=5000):
    """Largest erasure probability for which density evolution goes to zero."""
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if de_bec(lam, rho, mid, iters)[-1] < 1e-10:
            lo = mid
        else:
            hi = mid
    return lo


_gh_x, _gh_w = np.polynomial.hermite_e.hermegauss(80)      # Gauss-Hermite nodes


def phi_ga(m):
    """phi(m) = 1 - E[tanh(u/2)], u ~ N(m, 2m) (Chung et al., 2001), by quadrature."""
    m = np.atleast_1d(np.asarray(m, float))
    u = m[:, None] + np.sqrt(2 * m[:, None]) * _gh_x[None, :]
    val = 1 - (np.tanh(u / 2) * _gh_w).sum(axis=1) / np.sqrt(2 * np.pi)
    return np.where(m < 1e-10, 1.0, val)


_m_tab = np.logspace(-8, 3, 6000)                         # table for the inverse
_lphi_tab = np.log(np.maximum(phi_ga(_m_tab), 1e-300))


def phi_ga_inv(y):
    """Inverse of phi_ga by interpolation of ln(phi) on a log grid of m."""
    lm = np.interp(-np.log(max(y, 1e-300)), -_lphi_tab, np.log(_m_tab))
    return float(np.exp(lm))


def de_ga(lam, rho, sigma, iters=20000, big=60.0):
    """Gaussian-approximation density evolution on BI-AWGN. Returns True if the
    mean of the check-to-variable messages grows without bound."""
    m0, mu = 2 / sigma ** 2, 0.0
    for _ in range(iters):
        s = sum(c * phi_ga(m0 + (d - 1) * mu)[0] for d, c in lam.items())
        mu = sum(c * phi_ga_inv(1 - (1 - s) ** (d - 1)) for d, c in rho.items())
        if mu > big:
            return True
    return False


def threshold_ga(lam, rho):
    lo, hi = 0.3, 1.5
    for _ in range(30):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if de_ga(lam, rho, mid) else (lo, mid)
    return lo


def de_mc(dv, dc, sigma, iters=300, pop=100000, seed=1):
    """Monte-Carlo (population) density evolution for a (dv, dc)-regular code
    with sum-product decoding. Returns the bit error probability per iteration."""
    rng = np.random.default_rng(seed)
    ch = lambda: 2 / sigma ** 2 * (1 + sigma * rng.standard_normal(pop))
    v, pe = ch(), []
    for _ in range(iters):
        t = np.prod(np.tanh(v[rng.integers(0, pop, (dc - 1, pop))] / 2), axis=0)
        u = 2 * np.arctanh(np.clip(t, -1 + 1e-15, 1 - 1e-15))
        v = ch() + u[rng.integers(0, pop, (dv - 1, pop))].sum(axis=0)
        pe.append(np.mean(v < 0))
        if np.median(v) > 60:
            break
    return pe


def J(sigma):
    """Mutual information of a consistent Gaussian LLR with std sigma (EXIT charts)."""
    s = np.atleast_1d(np.asarray(sigma, float))
    m = s ** 2 / 2
    u = m[:, None] + s[:, None] * _gh_x[None, :]
    val = 1 - (np.logaddexp(0, -u) / np.log(2) * _gh_w).sum(axis=1) / np.sqrt(2 * np.pi)
    return np.where(s < 1e-8, 0.0, val)


def J_inv(I):
    lo, hi = 0.0, 60.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if J(mid)[0] < I else (lo, mid)
    return (lo + hi) / 2
# ==== [/de]


# ==== [sim]
# ---------------------------------------------------------------- 11. Monte-Carlo simulation
def simulate(encode, decode, k, R, ebn0_list, info=None, max_frames=100000,
             min_errors=100, batch=1000, seed=2026, stop_fer=0.0):
    """BER/FER versus Eb/N0. encode(u) -> (codeword, sent bits);
    decode(llr of sent bits) -> (decoded codeword, iterations); info = positions
    of the message bits inside the decoded codeword. The sweep stops after the
    first point with FER < stop_fer."""
    rng = np.random.default_rng(seed)
    res = {"ebn0": list(ebn0_list), "ber": [], "fer": [], "frames": [], "avg_iter": []}
    for eb in ebn0_list:
        bit_err = fr_err = frames = it_sum = 0
        while frames < max_frames and fr_err < min_errors:
            u = rng.integers(0, 2, (batch, k), dtype=np.uint8)
            c, x = encode(u)
            llr = awgn_llr(x, eb, R, rng)
            c_hat, iters = decode(llr)
            u_hat = c_hat[:, info] if info is not None else c_hat[:, :k]
            e = (u_hat != u).sum(axis=1)
            bit_err += e.sum(); fr_err += (e > 0).sum(); frames += batch; it_sum += iters.sum()
        res["ber"].append(bit_err / (frames * k)); res["fer"].append(fr_err / frames)
        res["frames"].append(frames); res["avg_iter"].append(it_sum / frames)
        if fr_err / frames < stop_fer:                     # low enough: stop the sweep
            res["ebn0"] = res["ebn0"][:len(res["fer"])]
            break
    return res
# ==== [/sim]
