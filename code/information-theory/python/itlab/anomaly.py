"""anomaly.py -- entropy- and KL-based anomaly detection (Lab 11).

A 'normal' user is modelled by a Markov chain over 4 actions
  0 = login, 1 = read, 2 = write, 3 = logout.
A baseline model is learned from training data. Each new window of W actions
is scored by the KL divergence D(empirical || baseline) of its action
frequencies and by the average log-loss under the Markov baseline.
A window is flagged when its score exceeds a threshold chosen from normal
validation windows (e.g. the 99th percentile).
"""
import random
from math import log2
from markov import estimate, simulate, stationary

ACTIONS = ["login", "read", "write", "logout"]
P_NORMAL = [[0.05, 0.80, 0.10, 0.05],
            [0.02, 0.60, 0.25, 0.13],
            [0.02, 0.50, 0.30, 0.18],
            [0.90, 0.04, 0.03, 0.03]]
P_ATTACK = [[0.60, 0.10, 0.05, 0.25],     # repeated login attempts
            [0.05, 0.30, 0.55, 0.10],     # unusual amount of writing
            [0.05, 0.20, 0.65, 0.10],
            [0.85, 0.05, 0.05, 0.05]]


def freq(window, S=4, alpha=0.5):
    c = [alpha] * S                         # add-alpha smoothing avoids log(0)
    for a in window:
        c[a] += 1
    t = sum(c)
    return [x / t for x in c]


def kl(p, q):
    return sum(pi * log2(pi / qi) for pi, qi in zip(p, q) if pi > 0)


def log_loss(window, P):
    """Average -log2 P(x_k | x_{k-1}) in bits per action (cross-entropy rate)."""
    return -sum(log2(P[a][b]) for a, b in zip(window, window[1:])) / (len(window) - 1)


def windows(seq, W):
    return [seq[i:i + W] for i in range(0, len(seq) - W + 1, W)]


def percentile(values, q):
    v = sorted(values)
    return v[min(len(v) - 1, int(q * len(v)))]


def run(W=50, seed=7):
    train = simulate(P_NORMAL, 20000, seed=seed)
    P_hat = estimate(train, 4, alpha=1.0)           # smoothed MLE baseline
    base = stationary(P_hat)
    valid = windows(simulate(P_NORMAL, 50 * W * 4, seed=seed + 1), W)
    test_norm = windows(simulate(P_NORMAL, 200 * W, seed=seed + 2), W)
    test_att = windows(simulate(P_ATTACK, 200 * W, seed=seed + 3), W)
    res = {}
    for name, score in [("KL", lambda w: kl(freq(w), base)),
                        ("log-loss", lambda w: log_loss(w, P_hat))]:
        thr = percentile([score(w) for w in valid], 0.99)
        fp = sum(score(w) > thr for w in test_norm) / len(test_norm)
        fn = sum(score(w) <= thr for w in test_att) / len(test_att)
        res[name] = (thr, fp, fn)
    return P_hat, base, res


if __name__ == "__main__":
    P_hat, base, res = run()
    print("baseline action frequencies:", [round(x, 3) for x in base])
    for name, (thr, fp, fn) in res.items():
        print("%-8s threshold %.3f  false positive rate %.3f  false negative rate %.3f"
              % (name, thr, fp, fn))
