"""Linear-programming design of irregular LDPC degree distributions for the BEC."""
import numpy as np
from scipy.optimize import linprog
from ldpc import threshold_bec, poly

def feasible(eps, dc, dmax, rate, grid=np.linspace(1e-4, 1, 600)):
    """Is there lambda (degrees 2..dmax) with design rate >= rate and
    eps*lambda(1-rho(1-x)) < x on the grid? rho(x) = x^(dc-1)."""
    degs = np.arange(2, dmax + 1)
    y = 1 - (1 - grid) ** (dc - 1)                       # 1 - rho(1 - x)
    A = eps * y[:, None] ** (degs[None, :] - 1)          # eps * y^(i-1)
    b = grid * (1 - 1e-6)
    A_ub = np.vstack([A, -1.0 / degs[None, :]])          # rate: sum lam_i/i >= (1/dc)/(1-rate)
    b_ub = np.r_[b, -(1 / dc) / (1 - rate)]
    A_ub = np.vstack([A_ub, np.r_[eps * (dc - 1), np.zeros(len(degs) - 1)]])   # stability
    b_ub = np.r_[b_ub, 1 - 1e-6]
    r = linprog(-1.0 / degs, A_ub=A_ub, b_ub=b_ub, A_eq=np.ones((1, len(degs))), b_eq=[1],
                bounds=[(0, 1)] * len(degs), method="highs")
    return r.status == 0, (dict(zip(degs.tolist(), r.x)) if r.status == 0 else None)

if __name__ == "__main__":
    for dc, dmax in [(6, 3), (6, 6), (7, 10), (8, 16), (9, 30)]:
        lo, hi, best = 0.3, 0.5, None
        for _ in range(30):
            mid = (lo + hi) / 2
            ok, lam = feasible(mid, dc, dmax, 0.5)
            if ok:
                lo, best = mid, lam
            else:
                hi = mid
        lam = {d: c for d, c in best.items() if c > 1e-4}
        s = sum(lam.values()); lam = {d: c / s for d, c in lam.items()}
        rate = 1 - (1 / dc) / sum(c / d for d, c in lam.items())
        th = threshold_bec(lam, {dc: 1.0})
        print(f"dc={dc} dmax={dmax}: LP eps={lo:.4f}  checked threshold={th:.4f}  rate={rate:.4f}  lambda=",
              {d: round(c, 4) for d, c in lam.items()})
