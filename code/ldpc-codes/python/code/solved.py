"""Numbers for the fully solved problems of Chapter 14."""
import numpy as np
from ldpc import *
np.set_printoptions(precision=4, suppress=True)
H6 = np.array([[1,1,0,1,0,0],[0,1,1,0,1,0],[1,0,0,0,1,1],[0,0,1,1,0,1]], dtype=np.uint8)
c6 = LDPCCode.from_H(H6)
print("P1:", 600/1200, 0.5*10**0.15, sigma_from_ebn0(1.5, 0.5), 2*(-0.4)/sigma_from_ebn0(1.5,0.5)**2)
print("P2:", (H6 @ np.array([1,1,1,0,0,1])) % 2)
H8 = np.array([[1,1,1,1,0,0,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0],[0,1,0,1,0,1,0,1]])
c8 = LDPCCode.from_H(H8)
print("P3: dv", c8.dv, "dc", c8.dc, "rank", gf2_rank(H8), "girth", girth(c8), "4cyc", count_4cycles(c8))
r = np.array([1,1,0,1,1,0]); s = (H6 @ r) % 2; print("P5: r", r, "s", s, "fails", H6.T @ s)
print("   BF:", bit_flip(c6, r, 10))
y = np.array([1,-1,-1,0,-1,1]); print("P6 peel:", peel(c6, y))
v = np.array([1.5,-2.0,0.8,3.0])
for i in range(4):
    o = np.delete(v, i); print(f"P7 edge {i}: exact {2*np.arctanh(np.prod(np.tanh(o/2))):.4f}  ms {np.prod(np.sign(o))*np.abs(o).min():.4f}")
Lch, inc = 0.6, np.array([-1.1, 0.4, 2.0]); app = Lch + inc.sum(); print("P8 app", app, "msgs", app - inc)
L = np.array([-0.5, 2.5, -1.5, 1.0, -2.0, 1.5])
print("P9 hard", (L<0).astype(int), "syndrome", c6.syndrome((L<0).astype(np.uint8)))
c2v = check_update(c6, L[c6.vn][None,:], "spa"); print("P9 c2v", c2v[0]); app = L + c6.var_sum(c2v)[0]
print("P9 app", app, "dec", (app<0).astype(int), "syn", c6.syndrome((app<0).astype(np.uint8)))
x, it, a = bp_decode(c6, L, 20, "spa"); print("P9 full", x, it, a)
c2v = check_update(c6, L[c6.vn][None,:], "nms"); print("P10 c2v", c2v[0]); app = L + c6.var_sum(c2v)[0]
print("P10 app", app, (app<0).astype(int)); print("P10 full", bp_decode(c6, L, 20, "nms"))
base = np.array([[0,2,4],[1,3,5]]); print("P12 girth", girth(qc_lift(base,6)), count_4cycles(qc_lift(base,6)))
for (i,j) in [(0,1),(0,2),(1,2)]:
    print("   cols", i, j, (base[0,i]-base[1,i]+base[1,j]-base[0,j]) % 6)
nr = NRCode(200, 600); print("P13", nr.bg, nr.Z, nr.ils, nr.Kp, nr.Kp-200, nr.code.m//nr.Z, nr.code.n//nr.Z)
print("P16", de_bec({4:1.0},{8:1.0},0.35,3)[:3], "final", de_bec({4:1.0},{8:1.0},0.35)[-1])
s = 0.9; m0 = 2/s**2; ph = phi_ga(m0)[0]; mu = phi_ga_inv(1-(1-ph)**5)
print("P19 m0", m0, "phi(m0)", ph, "mu1", mu, "mv2", m0 + 2*mu)
print("P20", 1-0.4*0.5**2, 0.9**5)
print("P25 BG", nr_select(6696, 0.5), 22*320-6696)
