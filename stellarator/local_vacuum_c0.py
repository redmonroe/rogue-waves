"""Vacuum test: why some runs with u . B = 0 end with an asymmetry of 1e-7 instead of 1e-13.

Ladder one seed from order 6 to order 9 (as local_vacuum.py --ladder), then study the end point
with the residual evaluated at 50 digits:
  1. refinement that only accepts steps reducing the residual: it sticks at about 7e-16;
  2. scans along the two directions in which the Jacobian nearly vanishes: the residual grows
     like the square of the move, and the stuck residual lies along one of them;
  3. plain Newton steps including those directions: the asymmetry halves at every step and the
     residual falls by four, the signature of a double root. It stops near 5e-10 only because
     the double-precision Jacobian cannot resolve smaller singular values.
Conclusion: the exact solution is symmetric; 1e-7 is the square root of the rounding level.

Usage:  OMP_NUM_THREADS=3 python local_vacuum_c0.py [seed]     (seed 13 by default; about 10 minutes)
"""
import sys

import mpmath as mp
import numpy as np
from local_vacuum import Vac, pad
seed, K0 = (int(sys.argv[1]) if len(sys.argv) > 1 else 13), 6
V = Vac(K0)
rng = np.random.default_rng(1000 * K0 + seed)
x, nr, it, ok = V.solve(rng.normal(size=V.nvar) * 0.4 ** V.var_deg, 300)
for K in range(7, 10):
    W = Vac(K); x, nr, it, ok = W.solve(pad(V, x, W), 300); V = W
print("K=9 residual", nr, "sym_B", V.sym_B(x), "C", x[-1])

mp.mp.dps = 50
Vd, Vh = V, Vac(9)
J = Vh.J
def mul(a, b):
    out = np.zeros(J.N, dtype=object)
    np.add.at(out, J.O, a[J.I] * b[J.J])
    return out
J.mul = mul
xh = np.array([mp.mpf(float(v)) for v in x], dtype=object)
print('refinement that only accepts steps which reduce the residual')
norm = lambda v: mp.sqrt(mp.fsum(t * t for t in v))
tof = lambda v: np.array([float(t) for t in v])
r = Vh.residual(xh); nr = norm(r)
for it in range(60):
    xf = tof(xh)
    U, s, Vt = np.linalg.svd(Vd.jacobian(xf), full_matrices=False)
    proj = U.T @ tof(r)
    best = None
    for cut in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8):
        keep = s > cut * s[0]
        for t in (1.0, 2.0):          # t = 2 is the right length at a double root
            xn = xh - t * (Vt[keep].T @ (proj[keep] / s[keep]))
            rn = Vh.residual(xn); nn = norm(rn)
            if best is None or nn < best[0]:
                best = (nn, xn, rn, cut, t)
    if best[0] >= nr:
        print("no further decrease"); break
    nr, xh, r = best[0], best[1], best[2]
    print(f"pass {it:2d}: residual {mp.nstr(nr, 3):>9s}  sym_B {Vd.sym_B(tof(xh)):.2e}  strain {Vd.strain(tof(xh)):.2e} "
          f" C {float(xh[-1]):+.1e}  cut {best[3]:.0e} t {best[4]:.0f}  small singular values: "
          + " ".join(f"{v:.0e}" for v in (s / s[0])[-33:-29]), flush=True)

# The residual is stuck. Scan along the two nearly-null directions of the Jacobian.
xf = tof(xh)
U, s, Vt = np.linalg.svd(Vd.jacobian(xf), full_matrices=False)
order = np.argsort(s)
near = [i for i in order if 1e-13 < s[i] / s[0] < 1e-8]
print("nearly-null singular values:", [f"{s[i]:.1e}" for i in near])
for i in near:
    v = Vt[i]
    for a in (-3e-6, -1e-6, -3e-7, -1e-7, 0, 1e-7, 3e-7, 1e-6, 3e-6):
        xn = xh + mp.mpf(a) * v
        print(f"  direction {i}: move {a:+.0e}  residual {mp.nstr(norm(Vh.residual(xn)), 3):>9s}  sym_B {Vd.sym_B(tof(xn)):.2e}")

print("fine scan")
r0 = Vh.residual(xh)
for i in near:
    v = Vt[i]
    q = (Vh.residual(xh + mp.mpf(1e-6) * v) + Vh.residual(xh - mp.mpf(1e-6) * v) - 2 * r0) / (2 * mp.mpf(1e-6) ** 2)
    lin = (Vh.residual(xh + mp.mpf(1e-6) * v) - Vh.residual(xh - mp.mpf(1e-6) * v)) / (2 * mp.mpf(1e-6))
    cosang = mp.fsum(a * b for a, b in zip(q, r0)) / (norm(q) * norm(r0))
    print(f"  direction {i}: |quadratic coefficient| {mp.nstr(norm(q), 3)}  |linear coefficient| {mp.nstr(norm(lin), 3)}"
          f"  cosine(quadratic term, residual) {mp.nstr(cosang, 3)}")
    for a in (-3e-8, -1e-8, -3e-9, 3e-9, 1e-8, 3e-8):
        print(f"     move {a:+.0e}  residual {mp.nstr(norm(Vh.residual(xh + mp.mpf(a) * v)), 4)}")
# how much of the stuck residual is outside the range of the Jacobian (with the two nearly-null directions dropped)?
keep = s > 1e-9 * s[0]
rf = tof(r0)
inside = U[:, keep] @ (U[:, keep].T @ rf)
print("stuck residual:", np.linalg.norm(rf), " part inside the range of the Jacobian:", np.linalg.norm(inside),
      " along the two nearly-null left vectors:", [float(abs(U[:, i] @ rf)) for i in near])

print("plain Newton steps including the nearly-null directions (no requirement that the residual decreases)")
for it in range(40):
    xf = tof(xh)
    U, s, Vt = np.linalg.svd(Vd.jacobian(xf), full_matrices=False)
    keep = s > 1e-14 * s[0]
    r = Vh.residual(xh)
    xh = xh - Vt[keep].T @ ((U.T @ tof(r))[keep] / s[keep])
    if it % 3 == 2:
        print(f"  step {it:2d}: residual {mp.nstr(norm(Vh.residual(xh)), 3):>9s}  sym_B {Vd.sym_B(tof(xh)):.2e}  "
              f"strain {Vd.strain(tof(xh)):.2e}  C {float(xh[-1]):+.1e}  nearly-null singular values "
              + " ".join(f"{v:.0e}" for v in s[-31:-29]), flush=True)
