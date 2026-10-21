"""markov.py -- Markov sources, entropy rate, estimation (Lab 9).

A first-order chain on states 0..S-1 has transition matrix P with
P[i][j] = Pr(next = j | current = i); every row sums to 1.
"""
import random
from math import log2


def stationary(P, iters=10000, tol=1e-14):
    """Stationary distribution pi = pi P by power iteration."""
    S = len(P)
    pi = [1.0 / S] * S
    for _ in range(iters):
        new = [sum(pi[i] * P[i][j] for i in range(S)) for j in range(S)]
        if max(abs(a - b) for a, b in zip(new, pi)) < tol:
            return new
        pi = new
    return pi


def row_entropy(row):
    return -sum(p * log2(p) for p in row if p > 0)


def entropy_rate(P, pi=None):
    """H = sum_i pi_i H(P[i]) bits per symbol (stationary first-order chain)."""
    if pi is None:
        pi = stationary(P)
    return sum(pi[i] * row_entropy(P[i]) for i in range(len(P)))


def sequence_probability(seq, P, p0):
    """P(x1..xn) = p0(x1) * prod P[x_{k-1}][x_k]."""
    pr = p0[seq[0]]
    for a, b in zip(seq, seq[1:]):
        pr *= P[a][b]
    return pr


def estimate(seq, S, alpha=0.0):
    """Maximum-likelihood transition matrix from counts N_ij / N_i
    (alpha > 0 adds Laplace smoothing)."""
    N = [[alpha] * S for _ in range(S)]
    for a, b in zip(seq, seq[1:]):
        N[a][b] += 1
    P = []
    for row in N:
        t = sum(row)
        P.append([x / t if t else 1.0 / S for x in row])
    return P


def simulate(P, n, start=0, seed=None):
    rng = random.Random(seed)
    s, out = start, [start]
    for _ in range(n - 1):
        u, c = rng.random(), 0.0
        for j, p in enumerate(P[s]):
            c += p
            if u < c:
                s = j
                break
        out.append(s)
    return out


if __name__ == "__main__":
    P = [[0.9, 0.1], [0.5, 0.5]]
    pi = stationary(P)
    print("stationary distribution:", [round(x, 4) for x in pi])
    print("entropy rate: %.4f bits/symbol" % entropy_rate(P))
    print("entropy of pi (memoryless view): %.4f" % row_entropy(pi))
    seq = simulate(P, 20000, seed=1)
    print("estimated P from 20000 samples:", [[round(x, 3) for x in r] for r in estimate(seq, 2)])
