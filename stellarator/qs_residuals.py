"""Quasisymmetry residuals for published closed-form equilibria.

For each field, at random points with 40-digit arithmetic:
  mhs    : max of |div B|, |(curl B) x B - grad p|, |B . grad psi|   (should be ~0: checks the transcription)
  triple : grad(psi) x grad(B^2) . grad(B . grad B^2), normalised      (0 iff weakly quasisymmetric with this psi)
  frob   : N . curl N with N = B x (grad B^2 x grad(B . grad B^2)), normalised
           (must vanish if ANY flux label makes the field weakly quasisymmetric,
            so it does not depend on which psi is used)
The axisymmetric member of each family is included as a control and must give zeros.
"""
import random
import sys
import time

import mpmath as mp
import sympy as sp

from fields import X, landreman_iota2, landreman_sheared

mp.mp.dps = 40


def grad(f):
    return sp.Matrix([sp.diff(f, v) for v in X])


def div(V):
    return sum(sp.diff(V[i], X[i]) for i in range(3))


def curl(V):
    x, y, z = X
    return sp.Matrix([
        sp.diff(V[2], y) - sp.diff(V[1], z),
        sp.diff(V[0], z) - sp.diff(V[2], x),
        sp.diff(V[1], x) - sp.diff(V[0], y),
    ])


def norm(v):
    return mp.sqrt(sum(abs(c) ** 2 for c in v))


def residuals(B, psi, p, in_domain, box, n_points=40, seed=0):
    B2 = (B.T * B)[0]
    BdB2 = (B.T * grad(B2))[0]
    gpsi, gB2, gBdB2 = grad(psi), grad(B2), grad(BdB2)
    N = B.cross(gB2.cross(gBdB2))
    mhs_exprs = [div(B)] + list(curl(B).cross(B) - grad(p)) + [(B.T * gpsi)[0]]
    f = sp.lambdify(X, [mhs_exprs, list(gpsi), list(gB2), list(gBdB2), list(N), list(B)], "mpmath", cse=True)
    N_fn = sp.lambdify(X, list(N), "mpmath", cse=True)

    def curl_N(pt):
        d = lambda i, j: mp.diff(lambda t: N_fn(*[t if k == j else pt[k] for k in range(3)])[i], pt[j])
        return [d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)]

    rng = random.Random(seed)
    worst = dict(mhs=mp.mpf(0), triple=mp.mpf(0), frob=mp.mpf(0))
    n = 0
    while n < n_points:
        pt = [mp.mpf(rng.uniform(lo, hi)) for lo, hi in box]
        try:
            if not in_domain(pt):
                continue
            mhs, g1, g2, g3, Nv, Bv = f(*pt)
        except (ZeroDivisionError, ValueError):
            continue
        Bn = norm(Bv)
        if Bn < 0.05:
            continue
        n += 1
        worst["mhs"] = max(worst["mhs"], max(abs(c) for c in mhs))
        triple = mp.re(g1[0] * (g2[1] * g3[2] - g2[2] * g3[1]) + g1[1] * (g2[2] * g3[0] - g2[0] * g3[2])
                       + g1[2] * (g2[0] * g3[1] - g2[1] * g3[0]))
        scale = norm(g1) * norm(g2) * norm(g3)
        if scale > 0:
            worst["triple"] = max(worst["triple"], abs(triple) / scale)
        cN = curl_N(pt)
        nN, ncN = norm(Nv), norm(cN)
        if nN > 0 and ncN > 0:
            worst["frob"] = max(worst["frob"], abs(mp.re(sum(a * b for a, b in zip(Nv, cN)))) / (nN * ncN))
    return worst


R = sp.Rational
CASES = [
    ("Landreman iota=2, eps=0 (axisymmetric control)", landreman_iota2(sp.Integer(1), sp.Integer(1)), [(-1.3, 1.3), (-1.3, 1.3), (-0.45, 0.45)]),
    ("Landreman iota=2, eps=1/4", landreman_iota2(sp.sqrt(R(5, 4)), sp.sqrt(R(3, 4))), [(-1.3, 1.3), (-1.3, 1.3), (-0.45, 0.45)]),
    ("Landreman iota=2, eps=1/2", landreman_iota2(sp.sqrt(R(3, 2)), sp.sqrt(R(1, 2))), [(-1.3, 1.3), (-1.3, 1.3), (-0.45, 0.45)]),
    ("Landreman sheared, eps=0, lam=1, S=1 (axisymmetric control)", landreman_sheared(sp.Integer(0), sp.Integer(1), sp.Integer(1)), [(-1.5, 1.5), (-1.5, 1.5), (-1.0, 1.0)]),
    ("Landreman sheared, eps=1/2, lam=1, S=1", landreman_sheared(R(1, 2), sp.Integer(1), sp.Integer(1)), [(-1.5, 1.5), (-1.5, 1.5), (-1.0, 1.0)]),
    ("Landreman sheared, eps=1, lam=2, S=1/2", landreman_sheared(sp.Integer(1), sp.Integer(2), R(1, 2)), [(-1.5, 1.5), (-1.5, 1.5), (-1.0, 1.0)]),
]

if __name__ == "__main__":
    n_points = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    print(f"{'field':62s} {'mhs':>10s} {'triple':>10s} {'frob':>10s}   ({n_points} points each)")
    for name, (B, psi, p, in_domain), box in CASES:
        t0 = time.time()
        w = residuals(B, psi, p, in_domain, box, n_points=n_points)
        print(f"{name:62s} {float(w['mhs']):10.1e} {float(w['triple']):10.1e} {float(w['frob']):10.1e}   [{time.time() - t0:.0f}s]", flush=True)
