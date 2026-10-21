"""entropy.py -- information measures (Labs 1 and 3/4 of the book).

All logarithms are base 2, so every result is in bits.
Probabilities equal to 0 are skipped, using the convention 0*log(0) = 0.
"""
from collections import Counter
from math import log2


def self_information(p):
    """I(x) = -log2 p(x): the surprise of an outcome with probability p."""
    if not 0 < p <= 1:
        raise ValueError("p must be in (0, 1]")
    return -log2(p)


def entropy(probs):
    """H(X) = -sum p log2 p for a list (or dict values) of probabilities."""
    if isinstance(probs, dict):
        probs = list(probs.values())
    if abs(sum(probs) - 1) > 1e-9:
        raise ValueError("probabilities must sum to 1")
    return -sum(p * log2(p) for p in probs if p > 0)


def binary_entropy(p):
    """H2(p) = -p log2 p - (1-p) log2 (1-p)."""
    return entropy([p, 1 - p])


def probs_from_text(text):
    """Relative frequency of every symbol, sorted by symbol."""
    counts = Counter(text)
    n = len(text)
    return {s: counts[s] / n for s in sorted(counts)}


# ---------------------------------------------------------------- two variables
def marginals(P):
    """P is a joint table P[x][y] (list of rows). Returns (p_x, p_y)."""
    px = [sum(row) for row in P]
    py = [sum(P[i][j] for i in range(len(P))) for j in range(len(P[0]))]
    return px, py


def joint_entropy(P):
    """H(X,Y) from a joint table."""
    return entropy([p for row in P for p in row])


def conditional_entropy(P):
    """H(Y|X) = H(X,Y) - H(X)  (rows = X, columns = Y)."""
    px, _ = marginals(P)
    return joint_entropy(P) - entropy(px)


def mutual_information(P):
    """I(X;Y) = H(X) + H(Y) - H(X,Y)."""
    px, py = marginals(P)
    return entropy(px) + entropy(py) - joint_entropy(P)


def kl_divergence(p, q):
    """D(p||q) = sum p log2(p/q). Infinite if q = 0 where p > 0."""
    d = 0.0
    for pi, qi in zip(p, q):
        if pi > 0:
            if qi == 0:
                return float("inf")
            d += pi * log2(pi / qi)
    return d


def min_entropy(probs):
    """H_inf(X) = -log2 max p: the best single guess of an attacker."""
    if isinstance(probs, dict):
        probs = list(probs.values())
    return -log2(max(probs))


if __name__ == "__main__":
    text = "ABRACADABRA"
    p = probs_from_text(text)
    print("text          :", text)
    print("probabilities :", {s: round(v, 4) for s, v in p.items()})
    print("H(X)          : %.4f bits/symbol" % entropy(p))
    print("max possible  : %.4f bits/symbol (log2 of %d symbols)" % (log2(len(p)), len(p)))
    print("H2(0.1)       : %.4f bits" % binary_entropy(0.1))
    P = [[0.25, 0.25], [0.0, 0.5]]            # joint table P[x][y]
    print("H(X,Y) = %.4f, H(Y|X) = %.4f, I(X;Y) = %.4f"
          % (joint_entropy(P), conditional_entropy(P), mutual_information(P)))
