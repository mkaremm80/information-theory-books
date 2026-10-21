"""channel.py -- discrete memoryless channels and capacity (Lab 10)."""
import random
from math import log2


def h2(p):
    return 0.0 if p in (0, 1) else -p * log2(p) - (1 - p) * log2(1 - p)


def bsc_capacity(p):
    return 1 - h2(p)


def bec_capacity(eps):
    return 1 - eps


def mutual_information(px, W):
    """I(X;Y) for input distribution px and channel matrix W[x][y] = P(y|x)."""
    ny = len(W[0])
    py = [sum(px[x] * W[x][y] for x in range(len(px))) for y in range(ny)]
    I = 0.0
    for x, pxv in enumerate(px):
        for y in range(ny):
            if pxv > 0 and W[x][y] > 0:
                I += pxv * W[x][y] * log2(W[x][y] / py[y])
    return I


def blahut_arimoto(W, iters=2000, tol=1e-12):
    """Capacity of any DMC: returns (C, best input distribution)."""
    nx, ny = len(W), len(W[0])
    px = [1.0 / nx] * nx
    for _ in range(iters):
        py = [sum(px[x] * W[x][y] for x in range(nx)) for y in range(ny)]
        D = [sum(W[x][y] * log2(W[x][y] / py[y]) for y in range(ny) if W[x][y] > 0)
             for x in range(nx)]
        new = [px[x] * 2 ** D[x] for x in range(nx)]
        s = sum(new)
        new = [v / s for v in new]
        if max(abs(a - b) for a, b in zip(new, px)) < tol:
            px = new
            break
        px = new
    return mutual_information(px, W), px


def bsc(bits, p, seed=None):
    """Flip each bit independently with probability p."""
    rng = random.Random(seed)
    return [b ^ (rng.random() < p) for b in bits]


if __name__ == "__main__":
    for p in (0.0, 0.01, 0.11, 0.5):
        print("BSC p = %.2f  C = %.4f bits/use" % (p, bsc_capacity(p)))
    W = [[0.9, 0.1], [0.2, 0.8]]           # a non-symmetric binary channel
    C, px = blahut_arimoto(W)
    print("Z-like channel W =", W, " C = %.4f, best P(X) = %s" % (C, [round(v, 4) for v in px]))
