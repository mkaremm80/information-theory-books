"""
polar.py  --  A small, readable polar-code library for the textbook
"Polar Codes: From Zero to Expert".

Conventions used everywhere in this book
----------------------------------------
* Kernel  F = [[1, 0],
               [1, 1]]
* Generator matrix  G_N = F^{(x)n}  (n-fold Kronecker power, NO bit reversal).
  This is the same convention used by 3GPP NR (TS 38.212).
* Encoding is a row-vector product over GF(2):   x = u G_N.
* Bit indices run from 0 to N-1 in natural order.
* BPSK mapping:  bit 0 -> +1,  bit 1 -> -1.
* LLR definition:  L = ln( P(bit = 0 | y) / P(bit = 1 | y) ).
  A positive LLR therefore favours bit 0.
"""
import numpy as np

# ---------------------------------------------------------------------------
# 1. Encoder
# ---------------------------------------------------------------------------
def generator_matrix(N):
    """Return G_N = F^{(x)n} as a 0/1 integer matrix."""
    F = np.array([[1, 0], [1, 1]], dtype=np.uint8)
    G = np.array([[1]], dtype=np.uint8)
    while G.shape[0] < N:
        G = np.kron(F, G) % 2          # F (x) G  -> [[G,0],[G,G]]
    return G


def polar_encode(u):
    """
    Butterfly (fast) encoder.  u may be a 1-D vector of length N or a 2-D
    array of shape (B, N) holding B blocks.  Returns x = u G_N (mod 2).

    Each of the log2(N) stages XORs the "upper" half of every butterfly
    into the "lower" half, exactly as in the butterfly diagram.
    """
    x = np.array(u, dtype=np.uint8, copy=True)
    one_d = (x.ndim == 1)
    if one_d:
        x = x[None, :]
    B, N = x.shape
    half = N // 2
    while half >= 1:                       # stages: span N/2, N/4, ..., 1
        x = x.reshape(B, -1, 2, half)      # groups of 2*half bits
        x[:, :, 0, :] ^= x[:, :, 1, :]     # upper part ^= lower part
        x = x.reshape(B, N)
        half //= 2
    return x[0] if one_d else x


# ---------------------------------------------------------------------------
# 2. Code construction (reliability of the bit-channels)
# ---------------------------------------------------------------------------
def bec_bhattacharyya(N, eps):
    """
    Exact Bhattacharyya parameters (= erasure probabilities) of the N
    synthetic bit-channels of a BEC(eps).
        Z(W-) = 2Z - Z^2      (worse  channel)
        Z(W+) = Z^2           (better channel)
    Natural index order, consistent with polar_encode().
    """
    if N == 1:
        return np.array([eps], dtype=float)
    zm = 2 * eps - eps ** 2
    zp = eps ** 2
    return np.concatenate([bec_bhattacharyya(N // 2, zm),
                           bec_bhattacharyya(N // 2, zp)])


# Table of log(phi(x)) for 0 < x < 10, computed by Gauss-Hermite quadrature of
#     phi(x) = 1 - E[ tanh(U/2) ],   U ~ N(x, 2x)
# (this is the exact definition used in Gaussian approximation).
_GH_T, _GH_W = np.polynomial.hermite.hermgauss(120)
_XGRID = np.logspace(-8, 1, 4000)
_LPHI_TAB = np.log(np.clip(1 - np.array(
    [np.sum(_GH_W * np.tanh((xx + 2 * np.sqrt(xx) * _GH_T) / 2)) / np.sqrt(np.pi)
     for xx in _XGRID]), 1e-300, None))


def _log_phi(x):
    """
    Natural log of the GA function phi(x).
      * 0 < x < 10 : exact value (table built by numerical integration)
      * x >= 10    : Chung's asymptotic form
                     phi(x) ~ sqrt(pi/x) exp(-x/4) (1 - 10/(7x))
    Working with log(phi) avoids underflow for very reliable channels.
    """
    if x <= _XGRID[0]:
        return np.log1p(-x / 2)            # phi(x) ~ 1 - x/2 near zero
    if x < 10:
        return float(np.interp(np.log(x), np.log(_XGRID), _LPHI_TAB))
    return 0.5 * np.log(np.pi / x) - x / 4 + np.log(1 - 10 / (7 * x))


def _ga_minus(m):
    """Mean LLR of the 'minus' (check-node) channel: phi^-1(1-(1-phi(m))^2)."""
    lp = _log_phi(m)
    p = np.exp(lp)
    target = lp + np.log(2 - p)          # log(2p - p^2) = log(1-(1-p)^2)
    lo, hi = 0.0, m                      # the minus channel is never better
    for _ in range(80):                  # bisection: log_phi is decreasing
        mid = 0.5 * (lo + hi)
        if _log_phi(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ga_llr_means_recursive(N, design_snr_db, rate):
    """
    Gaussian Approximation (GA) construction for BPSK over AWGN.
    Returns the mean LLR of every bit-channel (larger = more reliable),
    in natural index order (consistent with polar_encode()).
        minus child:  m- = phi^-1( 1 - (1 - phi(m))^2 )
        plus  child:  m+ = 2 m
    """
    esn0 = 10 ** (design_snr_db / 10) * rate
    m0 = 4 * esn0                        # channel LLR mean = 2/sigma^2
    def rec(n, m):
        if n == 1:
            return np.array([m])
        return np.concatenate([rec(n // 2, _ga_minus(m)), rec(n // 2, 2 * m)])
    return rec(N, m0)


def ga_construct(N, K, design_snr_db=2.0, rate=None):
    """Return a boolean mask: True = information position (K most reliable)."""
    if rate is None:
        rate = K / N
    m = ga_llr_means_recursive(N, design_snr_db, rate)
    info = np.zeros(N, dtype=bool)
    info[np.argsort(m)[-K:]] = True
    return info


def bec_construct(N, K, eps=0.5):
    z = bec_bhattacharyya(N, eps)
    info = np.zeros(N, dtype=bool)
    info[np.argsort(z)[:K]] = True
    return info


# ---------------------------------------------------------------------------
# 3. Channel: BPSK + AWGN
# ---------------------------------------------------------------------------
def bpsk(x):
    return 1.0 - 2.0 * x                       # 0 -> +1, 1 -> -1


def awgn(s, ebn0_db, rate, rng):
    """Add white Gaussian noise. Es = 1, Es/N0 = R * Eb/N0."""
    esn0 = 10 ** (ebn0_db / 10) * rate
    sigma = np.sqrt(1 / (2 * esn0))
    return s + sigma * rng.standard_normal(s.shape), sigma


def channel_llr(y, sigma):
    return 2 * y / sigma ** 2


# ---------------------------------------------------------------------------
# 4. f and g functions
# ---------------------------------------------------------------------------
def f_minsum(a, b):
    """f(a,b) ~= sign(a) sign(b) min(|a|,|b|)   (min-sum approximation)."""
    return np.sign(a) * np.sign(b) * np.minimum(np.abs(a), np.abs(b))


def f_exact(a, b):
    """f(a,b) = 2 atanh( tanh(a/2) tanh(b/2) )   (exact box-plus)."""
    t = np.tanh(a / 2) * np.tanh(b / 2)
    t = np.clip(t, -1 + 1e-15, 1 - 1e-15)
    return 2 * np.arctanh(t)


def g_func(a, b, c):
    """g(a,b,c) = b + (1-2c) a ,  c = partial-sum bit already decided."""
    return b + (1 - 2 * c.astype(float)) * a


# ---------------------------------------------------------------------------
# 5. SC / SCL / CA-SCL decoder (one recursive routine, vectorised)
# ---------------------------------------------------------------------------
class ListDecoder:
    """
    LLR-based successive-cancellation LIST decoder.
      * L = 1  gives the plain SC decoder.
      * crc  (optional) gives the CRC-aided SCL decoder.
    The decoder processes B received blocks in parallel (numpy arrays
    of shape (B, L, ...)) so that simulations run fast in pure Python.

    Path metric (min-sum form, Balatsoukas-Stimming et al. 2015):
        PM <- PM + |L|   if the bit decision disagrees with sign(L)
        PM <- PM         otherwise.     Smaller PM = more likely path.
    """

    def __init__(self, info_mask, L=1, crc=None, exact_f=False):
        self.info = np.asarray(info_mask, dtype=bool)
        self.N = self.info.size
        self.L = L
        self.crc = crc
        self.f = f_exact if exact_f else f_minsum

    # -- leaf: one bit decision for every path of every block --------------
    def _leaf(self, llr, idx):
        B, L = llr.shape
        pen0 = np.where(llr < 0, -llr, 0.0)       # cost of deciding 0
        pen1 = np.where(llr > 0, llr, 0.0)        # cost of deciding 1
        if not self.info[idx]:                    # frozen bit: always 0
            self.PM = self.PM + pen0
            u = np.zeros((B, L, 1), dtype=np.uint8)
            return u, u.copy(), np.tile(np.arange(L), (B, 1))
        # information bit: every path forks into two candidates
        cand = np.concatenate([self.PM + pen0, self.PM + pen1], axis=1)  # (B,2L)
        order = np.argsort(cand, axis=1, kind="stable")[:, :L]           # keep L best
        self.PM = np.take_along_axis(cand, order, axis=1)
        parent = order % L
        bit = (order // L).astype(np.uint8)[:, :, None]
        return bit, bit.copy(), parent

    # -- recursive node --------------------------------------------------------
    def _node(self, llr, start):
        """llr: (B, L, n).  Returns (codeword c, info vector u, parent map)."""
        n = llr.shape[2]
        if n == 1:
            return self._leaf(llr[:, :, 0], start)
        h = n // 2
        a, b = llr[:, :, :h], llr[:, :, h:]
        # left child: LLRs from f
        c1, u1, p1 = self._node(self.f(a, b), start)
        # right child: LLRs from g, using the surviving parents' LLRs
        a1 = np.take_along_axis(a, p1[:, :, None], axis=1)
        b1 = np.take_along_axis(b, p1[:, :, None], axis=1)
        c2, u2, p2 = self._node(g_func(a1, b1, c1), start + h)
        # re-order left results according to the right child's survivors
        c1 = np.take_along_axis(c1, p2[:, :, None], axis=1)
        u1 = np.take_along_axis(u1, p2[:, :, None], axis=1)
        parent = np.take_along_axis(p1, p2, axis=1)
        return np.concatenate([c1 ^ c2, c2], axis=2), \
               np.concatenate([u1, u2], axis=2), parent

    def decode(self, llr):
        """
        llr : (B, N) channel LLRs.
        Returns u_hat (B, N) chosen by: best PM (SC/SCL) or
        best PM among CRC-passing paths (CA-SCL).
        """
        llr = np.atleast_2d(llr)
        B = llr.shape[0]
        self.PM = np.full((B, self.L), np.inf)
        self.PM[:, 0] = 0.0                      # only one active path at start
        lam = np.repeat(llr[:, None, :], self.L, axis=1)
        _, u, _ = self._node(lam, 0)
        self.last_u_list, self.last_pm = u, self.PM
        # paths are already sorted by PM (best first) after the last split
        if self.crc is None:
            return u[:, 0, :]
        payload = u[:, :, self.info]             # (B, L, K)
        ok = self.crc.check(payload) & np.isfinite(self.PM)
        pick = np.where(ok.any(axis=1), ok.argmax(axis=1), 0)
        return u[np.arange(B), pick, :]


# ---------------------------------------------------------------------------
# 6. CRC (as a linear map over GF(2))
# ---------------------------------------------------------------------------
CRC_POLYS = {
    # degree-list form (exponents with coefficient 1), from 3GPP TS 38.212 5.1
    "CRC6":   [6, 5, 0],
    "CRC11":  [11, 10, 9, 5, 0],
    "CRC24C": [24, 23, 21, 20, 17, 15, 13, 12, 8, 4, 2, 1, 0],
}


def crc_bits(msg, exps):
    """Plain long-division CRC (register starts at zero). msg: 1-D 0/1."""
    r = exps[0]
    g = np.zeros(r + 1, dtype=np.uint8)
    for e in exps:
        g[r - e] = 1                        # g[0] is the x^r coefficient
    reg = np.concatenate([np.array(msg, dtype=np.uint8), np.zeros(r, np.uint8)])
    for i in range(len(msg)):
        if reg[i]:
            reg[i:i + r + 1] ^= g
    return reg[-r:]


class CRC:
    """CRC that appends r parity bits: [message | crc]."""

    def __init__(self, name, k_msg):
        self.exps = CRC_POLYS[name]
        self.r = self.exps[0]
        self.k = k_msg
        # Because the CRC is linear, crc(m) = m @ P  (mod 2). Build P once.
        P = np.zeros((k_msg, self.r), dtype=np.uint8)
        for i in range(k_msg):
            e = np.zeros(k_msg, dtype=np.uint8)
            e[i] = 1
            P[i] = crc_bits(e, self.exps)
        self.P = P

    def attach(self, m):
        m = np.atleast_2d(m).astype(np.uint8)
        return np.concatenate([m, (m.astype(int) @ self.P) % 2], axis=1).astype(np.uint8)

    def check(self, word):
        """word (..., k+r) -> True where the CRC is satisfied."""
        m, c = word[..., :self.k].astype(int), word[..., self.k:]
        return np.all((m @ self.P) % 2 == c, axis=-1)


# ---------------------------------------------------------------------------
# 7. Monte-Carlo simulation helper
# ---------------------------------------------------------------------------
def simulate(N, K, ebn0_list, L=1, crc_name=None, design_snr_db=None,
             max_frames=20000, min_errors=200, batch=500, seed=1):
    """
    Returns dict with BER and BLER (on the MESSAGE bits) for each Eb/N0.
    K counts all bits placed on information positions (message + CRC).
    The rate used for Eb/N0 is the message rate (K - r)/N.
    """
    rng = np.random.default_rng(seed)
    r = 0 if crc_name is None else CRC_POLYS[crc_name][0]
    k_msg = K - r
    rate = k_msg / N
    dsnr = design_snr_db if design_snr_db is not None else 2.0
    info = ga_construct(N, K, dsnr, rate=K / N)
    crc = CRC(crc_name, k_msg) if crc_name else None
    dec = ListDecoder(info, L=L, crc=crc)
    ber, bler = [], []
    for eb in ebn0_list:
        bit_err = blk_err = frames = 0
        while frames < max_frames and blk_err < min_errors:
            m = rng.integers(0, 2, size=(batch, k_msg), dtype=np.uint8)
            payload = crc.attach(m) if crc else m
            u = np.zeros((batch, N), dtype=np.uint8)
            u[:, info] = payload
            x = polar_encode(u)
            y, sigma = awgn(bpsk(x), eb, rate, rng)
            u_hat = dec.decode(channel_llr(y, sigma))
            m_hat = u_hat[:, info][:, :k_msg]
            e = (m_hat != m)
            bit_err += e.sum()
            blk_err += e.any(axis=1).sum()
            frames += batch
        ber.append(bit_err / (frames * k_msg))
        bler.append(blk_err / frames)
    return {"ebn0": list(ebn0_list), "ber": ber, "bler": bler}
