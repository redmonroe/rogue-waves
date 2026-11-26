"""Run the local Taylor-series test on a known exact non-symmetric solution.

The twisted elliptical pinch (planar_family.py, tan(theta/2) = y) is an exact quasisymmetric
equilibrium with no rigid symmetry. Its Taylor series about a generic point, moved into the
normalisation that local_jets.py uses, must

  1. make every residual of local_jets.py vanish (a test of that code against a non-symmetric
     exact solution; until now it had only been tested on a symmetric one), and
  2. keep a positive symmetry measure at every order, and
  3. keep free parameters at every degree when counted by local_jets_project.py (the family
     has a free function theta(y)), instead of collapsing like the branch followed there.

Usage:  python local_planar_check.py [Kmax]
"""
import sys
import time

import numpy as np

from local_jets_count import Variant
from local_jets_project import projected_counts

ALPHA, BETA = 1.0, 0.1
X0 = np.array([0.114, -0.473, 0.111])   # chosen so the series is well scaled (largest coefficient about 3)


def closed_form(r):
    """B, u, psi, p of the rational twisted pinch at a point r (plain floats)."""
    x, y, z = r
    W = 1 + y * y
    c, s = (1 - y * y) / W, 2 * y / W
    xp, zp = c * x + s * z, -s * x + c * z
    b0 = (-2 * BETA * zp, 2 * ALPHA * xp)
    B = np.array([c * b0[0] - s * b0[1], 0.0, s * b0[0] + c * b0[1]])
    u = np.array([-2 * z / W, 1.0, 2 * x / W])
    return B, u


def series(J, Q, scale):
    """Taylor series (as local_jets coefficient vectors) of the closed form in new coordinates xi,
    with r = X0 + scale * Q^T xi."""
    m = J.mul
    one = np.zeros(J.N)
    one[J.i0] = 1.0
    r = []
    for i in range(3):
        a = X0[i] * one
        for j in range(3):
            a[J.lin[j]] += scale * Q[j, i]
        r.append(a)
    x, y, z = r
    W = one + m(y, y)
    Winv = one / W[J.i0]
    for _ in range(8):  # Newton iteration for the reciprocal of a power series
        Winv = m(Winv, 2 * one - m(W, Winv))
    c, s = m(one - m(y, y), Winv), 2 * m(y, Winv)
    xp, zp = m(c, x) + m(s, z), m(c, z) - m(s, x)
    b0x, b0z = -2 * BETA * zp, 2 * ALPHA * xp
    B = [m(c, b0x) - m(s, b0z), np.zeros(J.N), m(s, b0x) + m(c, b0z)]
    u = [-2 * m(z, Winv), one.copy(), 2 * m(x, Winv)]
    psi = -(ALPHA * m(xp, xp) + BETA * m(zp, zp))
    p = 2 * (ALPHA + BETA) * psi
    return B, u, psi, p


def normalised(J):
    """The solution in the frame and units of local_jets.py:
    B(0) = (0,0,1), u(0) = (1,0,c), d_z B_3(0) = 1."""
    B0, u0 = closed_form(X0)
    nB = np.linalg.norm(B0)
    e3 = B0 / nB
    up = u0 - (u0 @ e3) * e3
    lam = 1.0 / np.linalg.norm(up)
    e1 = up * lam
    Q = np.array([e1, np.cross(e3, e1), e3])
    rot = lambda V: [sum(Q[i, j] * V[j] for j in range(3)) for i in range(3)]
    B, _, _, _ = series(J, Q, 1.0)
    scale = nB / rot(B)[2][J.lin[2]]          # makes d_z B_3(0) = 1 after dividing B by |B(0)|
    B, u, psi, p = series(J, Q, scale)
    fields = [b / nB for b in rot(B)] + [lam * v for v in rot(u)] + [lam / (nB * scale) * psi, p / nB ** 2]
    return np.concatenate([f[mk] for f, mk in zip(fields, J.masks)])


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    for K in range(3, Kmax + 1):
        t0 = time.time()
        J = Variant(K, "qs")
        x = normalised(J)
        res = np.linalg.norm(J.residual(x))
        dims, gap = projected_counts(J, x, min(9, K - 1))
        print(f"K={K:2d}  residual {res:7.1e}  sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  "
              f"freedom in B, k=1.. : " + " ".join(f"{v:4d}" for v in dims)
              + f"   min log10 gap {np.log10(gap):4.1f}   [{time.time() - t0:.0f}s]", flush=True)
