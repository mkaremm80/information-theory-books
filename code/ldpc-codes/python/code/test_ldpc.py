"""Unit tests of ldpc.py (run: python3 test_ldpc.py)."""
import itertools, numpy as np
from ldpc import *

H6 = np.array([[1, 1, 0, 1, 0, 0], [0, 1, 1, 0, 1, 0], [1, 0, 0, 0, 1, 1], [0, 0, 1, 1, 0, 1]], dtype=np.uint8)
rng = np.random.default_rng(1)
results = []

def test(name):
    def deco(f):
        f(); results.append(name); print("PASS", name)
    return deco

@test("GF(2): rank, inverse and generator matrix")
def _():
    assert gf2_rank(H6) == 3
    for _ in range(20):
        A = rng.integers(0, 2, (8, 8)).astype(np.uint8)
        if gf2_rank(A) == 8:
            assert ((A.astype(int) @ gf2_inv(A)) % 2 == np.eye(8)).all()
    for H in [H6, gallager_H(60, 3, 6, seed=2), peg(96, 48, 3, seed=0).H()]:
        G, info = generator_from_H(H)
        assert not ((H.astype(int) @ G.T) % 2).any()
        assert G.shape[0] == H.shape[1] - gf2_rank(H) == gf2_rank(G)
        assert (G[:, info] == np.eye(G.shape[0])).all()

@test("Code object: degrees, H and syndromes")
def _():
    c = LDPCCode.from_H(H6)
    assert (c.H() == H6).all() and c.E == 12 and (c.dv == 2).all() and (c.dc == 3).all()
    cws = [w for w in itertools.product([0, 1], repeat=6) if not ((H6 @ np.array(w)) % 2).any()]
    assert len(cws) == 8 and not c.syndrome(np.array(cws)).any()

@test("Girth and 4-cycles")
def _():
    assert girth(LDPCCode.from_H(H6)) == 6
    assert girth(qc_lift([[0, 0, 0, 0], [0, 1, 2, 3]], 5)) == 8
    q = qc_lift([[0, 0, 0, 0], [0, 1, 1, 3]], 5)
    assert girth(q) == 4 and count_4cycles(q) == 5
    tree = LDPCCode.from_H(np.array([[1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1]]))
    assert girth(tree) == np.inf

@test("5G NR tables, parameter selection and encoder")
def _():
    for ils in range(8):
        assert (nr_base_graph(1, ils) >= 0).sum() == 316 and (nr_base_graph(2, ils) >= 0).sum() == 197
    assert nr_select(120, 1 / 3)[:2] == (2, 20) and nr_select(1024, 0.5)[:2] == (2, 104)
    assert nr_select(1024, 8 / 9)[:2] == (1, 48) and nr_select(8448, 1 / 3)[:2] == (1, 384)
    for K, E in [(40, 200), (120, 360), (500, 1250), (1024, 2048), (1024, 1152), (3840, 19200), (8448, 25344), (100, 3000)]:
        nr = NRCode(K, E)
        c, x = nr.encode(rng.integers(0, 2, (3, K)))
        assert not nr.code.syndrome(c).any() and x.shape == (3, E)
        assert not c[:, nr.filler].any()

@test("Check-node update equals the brute-force formulas")
def _():
    c = LDPCCode.from_H(H6)
    v2c = rng.normal(0, 3, (50, c.E))
    spa = check_update(c, v2c, "spa"); ms = check_update(c, v2c, "ms")
    for e in range(c.E):
        others = [f for f in range(c.E) if c.cn[f] == c.cn[e] and f != e]
        t = np.prod(np.tanh(v2c[:, others] / 2), axis=1)
        assert np.allclose(spa[:, e], 2 * np.arctanh(t), atol=1e-8)
        assert np.allclose(ms[:, e], np.prod(np.sign(v2c[:, others]), 1) * np.abs(v2c[:, others]).min(1))
    assert np.allclose(check_update(c, v2c, "nms", alpha=0.75), 0.75 * ms)

@test("Sum-product is exact (bit-MAP) on a cycle-free graph")
def _():
    Ht = np.array([[1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1]], dtype=np.uint8)
    ct = LDPCCode.from_H(Ht)
    cw = np.array([w for w in itertools.product([0, 1], repeat=7) if not ((Ht @ np.array(w)) % 2).any()])
    for _ in range(100):
        L = rng.normal(1, 2, 7)
        app = bp_decode(ct, L, 10, "spa", early_stop=False)[2][0]
        ll = (np.where(cw == 0, 1, -1) * L / 2).sum(1)
        mp = [np.logaddexp.reduce(ll[cw[:, i] == 0]) - np.logaddexp.reduce(ll[cw[:, i] == 1]) for i in range(7)]
        assert np.abs(app - mp).max() < 1e-9

@test("Decoders: noiseless frames and high-SNR frames (flooding and layered)")
def _():
    code = peg(504, 252, 3, seed=3); G, info = generator_from_H(code.H()); make_layers(code)
    c = encode_G(G, rng.integers(0, 2, (200, G.shape[0]))).astype(np.uint8)
    L = 10.0 * bpsk(c)
    for rule in ["spa", "ms", "nms", "oms"]:
        x, it, _ = bp_decode(code, L, 5, rule); assert (x == c).all() and (it == 1).all()
    L = awgn_llr(c, 4.0, 0.5, np.random.default_rng(9))
    for f in [lambda: bp_decode(code, L, 50, "spa"), lambda: layered_decode(code, L, 50, "nms")]:
        assert (f()[0] == c).all()

@test("Peeling and bit flipping")
def _():
    code = peg(96, 48, 3, seed=0); G, _ = generator_from_H(code.H())
    c = encode_G(G, rng.integers(0, 2, (50, G.shape[0]))).astype(int)
    y = c.copy(); y[:, :3] = -1                                  # 3 erasures
    out, _ = peel(code, y)
    assert ((out == c) | (out < 0)).all()                        # never wrong
    assert girth(code) >= 6
    for i in range(96):                                          # every single error corrected
        r = c[:1].copy(); r[0, i] ^= 1
        assert (bit_flip(code, r, 10)[0] == c[:1]).all()

@test("Density evolution thresholds")
def _():
    assert abs(threshold_bec({3: 1.0}, {6: 1.0}) - 0.4294) < 1e-4
    assert abs(threshold_bec({3: 1.0}, {4: 1.0}) - 0.6474) < 1e-4
    assert 0.870 < threshold_ga({3: 1.0}, {6: 1.0}) < 0.874
    assert de_mc(3, 6, 0.86, iters=100, pop=50000)[-1] == 0.0

@test("Monte-Carlo simulation runs and decodes at high SNR")
def _():
    nr = NRCode(256, 512); make_layers(nr.code, nr.Z)
    dec = lambda L: layered_decode(nr.code, nr.llr_in(L), 20, "nms")[:2]
    r = simulate(nr.encode, dec, 256, 0.5, [4.0], info=np.arange(256), max_frames=400, batch=200)
    assert r["fer"][0] == 0.0 and r["frames"][0] == 400

print(f"\n{len(results)} tests passed")
