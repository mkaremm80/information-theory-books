"""test_itlab.py -- unit tests for every lab.  Run:  python test_itlab.py"""
import itertools, random, unittest
from fractions import Fraction
from math import log2
import entropy as E, shannon_fano as sf, huffman as hf, arithmetic as ac, lzw
import hamming as hm, repetition as rp, bch, markov as mk, channel as ch

rng = random.Random(1)


def random_text(n, alphabet):
    return "".join(rng.choice(alphabet) for _ in range(n))


def is_prefix_free(codes):
    cw = list(codes.values())
    return not any(a != b and b.startswith(a) for a in cw for b in cw)


class TestEntropy(unittest.TestCase):
    def test_values(self):
        self.assertAlmostEqual(E.entropy([0.5, 0.5]), 1.0)
        self.assertAlmostEqual(E.entropy([0.25] * 4), 2.0)
        self.assertAlmostEqual(E.entropy([1.0]), 0.0)
        self.assertAlmostEqual(E.binary_entropy(0.11), 0.4999, places=3)

    def test_identities(self):
        P = [[0.1, 0.2, 0.1], [0.3, 0.05, 0.25]]
        px, py = E.marginals(P)
        I = E.mutual_information(P)
        self.assertGreaterEqual(I, 0)
        self.assertAlmostEqual(E.joint_entropy(P), E.entropy(px) + E.conditional_entropy(P))
        # independence -> I = 0
        Q = [[a * b for b in [0.3, 0.7]] for a in [0.4, 0.6]]
        self.assertAlmostEqual(E.mutual_information(Q), 0.0)
        # I = D(p_xy || p_x p_y)
        flat = [p for r in P for p in r]
        prod = [a * b for a in px for b in py]
        self.assertAlmostEqual(I, E.kl_divergence(flat, prod))


class TestSourceCodes(unittest.TestCase):
    def test_prefix_kraft_bounds(self):
        for _ in range(200):
            k = rng.randint(2, 12)
            w = [rng.random() + 0.01 for _ in range(k)]
            probs = {chr(65 + i): x / sum(w) for i, x in enumerate(w)}
            H = E.entropy(probs)
            for coder in (sf.shannon_fano, hf.huffman):
                codes = coder(probs)
                self.assertTrue(is_prefix_free(codes))
                self.assertLessEqual(sum(2 ** -len(c) for c in codes.values()), 1 + 1e-12)
                L = sum(probs[s] * len(codes[s]) for s in probs)
                self.assertGreaterEqual(L, H - 1e-9)
                self.assertLess(L, H + 1)
            Lh = sum(probs[s] * len(hf.huffman(probs)[s]) for s in probs)
            Ls = sum(probs[s] * len(sf.shannon_fano(probs)[s]) for s in probs)
            self.assertLessEqual(Lh, Ls + 1e-12)            # Huffman is optimal

    def test_round_trip(self):
        for _ in range(50):
            t = random_text(rng.randint(1, 300), "ABCDEFG")
            p = E.probs_from_text(t)
            for m in (sf, hf):
                c = (sf.shannon_fano if m is sf else hf.huffman)(p)
                self.assertEqual(m.decode(m.encode(t, c), c), t)


class TestArithmetic(unittest.TestCase):
    def test_exact(self):
        probs = {"A": Fraction(1, 2), "B": Fraction(1, 4), "C": Fraction(1, 4)}
        for _ in range(100):
            t = random_text(rng.randint(1, 20), "ABC")
            lo, hi, bits = ac.encode_exact(t, probs)
            self.assertEqual(ac.decode_exact(bits, probs, len(t)), t)
            self.assertLessEqual(len(bits), -log2(hi - lo) + 2)

    def test_integer(self):
        from collections import Counter
        for _ in range(100):
            t = random_text(rng.randint(1, 2000), "ABCDE" if rng.random() < .5 else "AAAAAAB")
            c = Counter(t)
            self.assertEqual(ac.decode_int(ac.encode_int(t, c), c, len(t)), t)


class TestLZW(unittest.TestCase):
    def test_round_trip(self):
        for _ in range(200):
            t = random_text(rng.randint(1, 500), "AB" if rng.random() < .5 else "ABCD")
            self.assertEqual(lzw.decode(lzw.encode(t, "ABCD"), "ABCD"), t)
        self.assertEqual(lzw.encode("ABABABA", "AB"), [0, 1, 2, 4])   # special case
        t = "TOBEORNOTTOBEORTOBEORNOT"
        self.assertEqual(lzw.decode(lzw.encode(t)), t)


class TestChannelCodes(unittest.TestCase):
    def test_hamming(self):
        for d in itertools.product([0, 1], repeat=4):
            c = hm.encode74(list(d))
            Gc = [sum(d[i] * hm.G74[i][j] for i in range(4)) % 2 for j in range(7)]
            self.assertEqual(c, Gc)
            self.assertEqual(hm.syndrome74(c), 0)
            for e in range(7):                     # every single error corrected
                r = list(c); r[e] ^= 1
                self.assertEqual(hm.decode74(r)[0], list(d))
            for e1, e2 in itertools.combinations(range(8), 2):   # SECDED detects doubles
                r = hm.encode84(list(d)); r[e1] ^= 1; r[e2] ^= 1
                self.assertEqual(hm.decode84(r)[1], "double error detected")

    def test_hamming_two_errors_fail(self):
        c = hm.encode74([1, 0, 1, 1])
        r = list(c); r[1] ^= 1; r[5] ^= 1
        self.assertNotEqual(hm.decode74(r)[0], [1, 0, 1, 1])

    def test_repetition(self):
        for b in (0, 1):
            c = rp.encode([b])
            for k in range(6):
                for pos in itertools.combinations(range(5), k):
                    ok = rp.decode(rp.flip(c, pos))[0] == [b]
                    self.assertEqual(ok, k <= 2)

    def test_bch(self):
        for _ in range(30):
            m = [rng.randint(0, 1) for _ in range(7)]
            c = bch.encode(m)
            self.assertEqual(bch.syndromes(c), (0, 0))
            for k in (1, 2):
                for pos in itertools.combinations(range(15), k):
                    r = list(c)
                    for p in pos:
                        r[p] ^= 1
                    mh, chat, found = bch.decode(r)
                    self.assertEqual(chat, c)
                    self.assertEqual(sorted(found), list(pos))
        # minimum distance = 5: check all 128 codewords
        w = [sum(bch.encode(list(m))) for m in itertools.product([0, 1], repeat=7)]
        self.assertEqual(min(x for x in w if x > 0), 5)


class TestMarkovChannel(unittest.TestCase):
    def test_markov(self):
        P = [[0.9, 0.1], [0.5, 0.5]]
        pi = mk.stationary(P)
        self.assertAlmostEqual(pi[0], 5 / 6, places=9)
        self.assertAlmostEqual(mk.entropy_rate(P), 5 / 6 * E.binary_entropy(.1) + 1 / 6, places=9)

    def test_capacity(self):
        self.assertAlmostEqual(ch.bsc_capacity(0.11), 1 - E.binary_entropy(0.11))
        C, px = ch.blahut_arimoto([[0.89, 0.11], [0.11, 0.89]])
        self.assertAlmostEqual(C, ch.bsc_capacity(0.11), places=8)
        C, px = ch.blahut_arimoto([[1 - .3, .3, 0], [0, .3, .7]])     # BEC(0.3)
        self.assertAlmostEqual(C, 0.7, places=8)


if __name__ == "__main__":
    unittest.main(verbosity=1)
